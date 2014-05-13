using GaussDCA
using NLopt

immutable type PlmAlg
    verbose::Bool
    epsconv::Float64
    maxit::Int
end

immutable type PlmVar
    N::Int
    M::Int
    q::Int    
    q2::Int
    lambdaJ::Float64
    lambdaH::Float64
    Z::SharedArray{Int,2}
    W::SharedArray{Float64,1}

    function PlmVar(N,M,q,q2,lambdaJ, lambdaH, Z,W)
        sZ = SharedArray(Int,size(Z))
        sZ[:] = Z
        sW = SharedArray(Float64,size(W))
        sW[:] = W
        new(N,M,q,q2,lambdaJ, lambdaH, sZ,sW)
    end
end

function plmdca(filename::String; 
                max_gap_fraction::Real = 0.9, 
                theta = :auto, 
                lambdaJ::Real=0.005, 
                lambdaH::Real=0.01,
                epsconv::Real=1.0e-5,
                maxit::Int=1000,
                verbose::Bool=true)
    
    W,Z,N,M,q = ReadFasta(filename,max_gap_fraction, theta)

    plmalg = PlmAlg(verbose, epsconv ,maxit)
    plmvar = PlmVar(N,M,q,q*q,lambdaJ,lambdaH,Z,W)

    Jmat = MinimizePL(plmalg, plmvar)

    J, FN = ComputeScore(Jmat, plmvar)



    return J, FN
end
    

function ComputeScore(Jmat::Array{Float64,2}, var::PlmVar)

    q = var.q
    N = var.N

    JJ=reshape(Jmat[1:end-q,:], q,q,N-1,N)
    Jtemp1=zeros( q,q,int(N*(N-1)/2))
    Jtemp2=zeros( q,q,int(N*(N-1)/2))
    l = 1

    for i=1:(N-1)
        for j=(i+1):N
            Jtemp1[:,:,l]=JJ[:,:,j-1,i]; #J_ij as estimated from from g_i.
            Jtemp2[:,:,l]=JJ[:,:,i,j]'; #J_ij as estimated from from g_j.
            l=l+1;
        end
    end
    
    J1=zeros(q,q,int(N*(N-1)/2))
    J2=zeros(q,q,int(N*(N-1)/2))

    for l=1:int(N*(N-1)/2)
        J1[:,:,l] = Jtemp1[:,:,l]-repmat(mean(Jtemp1[:,:,l],1),q,1)-repmat(mean(Jtemp1[:,:,l],2),1,q) .+ mean(Jtemp1[:,:,l])
        J2[:,:,l] = Jtemp2[:,:,l]-repmat(mean(Jtemp2[:,:,l],1),q,1)-repmat(mean(Jtemp2[:,:,l],2),1,q) .+ mean(Jtemp2[:,:,l])
    end
    J = 0.5 * ( J1 + J2 )
    FN = zeros(Float64, N,N)
    l = 1
    for i=1:N-1
        for j=i+1:N
            FN[i,j] = vecnorm(J[:,:,l],2)
            FN[j,i] =FN[i,j]
            l+=1
        end
    end
    FN=GaussDCA.correct_APC(FN)
    return J,FN
end

function MinimizePL(alg::PlmAlg, var::PlmVar)

    x0 = zeros(Float64,(var.N - 1) * var.q2 + var.q)

    Jmat = @parallel hcat for site=1:var.N #1:12
        function f(x::Vector, g::Vector)
            if g === nothing
                g = zeros(Float64, length(x))
            end
            return PLsiteAndGrad!(x, site, g, var)            
        end
        
        opt = Opt(:LD_LBFGS, length(x0))
        ftol_abs!(opt, alg.epsconv)
        maxeval!(opt, alg.maxit)
        min_objective!(opt, f)
        elapstime = @elapsed  (minf, minx, ret) = optimize(opt, x0)
        if alg.verbose 
            @printf("site = %d\t pl = %.4f\t time = %.4f\n", site, minf, elapstime)
        end
        minx
    end 
    return Jmat
end

function ReadFasta(filename::String,max_gap_fraction::Real, theta)
    Z = GaussDCA.read_fasta_alignment(filename, max_gap_fraction)
    N, M = size(Z)
    q = int(maximum(Z))
    
    q > 32 && error("parameter q=$q is too big (max 31 is allowed)")
    _, _, Meff, W = GaussDCA.compute_new_frequencies(Z, theta)
    W  ./= Meff  
    Zint = int( Z )
    return W, Zint,N,M,q
end

function PLsiteAndGrad!(vecJ::Array{Float64,1}, site::Int, grad::Array{Float64,1}, plmvar::PlmVar)

    LL = length(vecJ)
    q2 = plmvar.q2
    q = plmvar.q
    N = plmvar.N
    M = plmvar.M
    Z = sdata(plmvar.Z)
    W = sdata(plmvar.W)

    for i=1:LL-q
        grad[i] = 2.0 * plmvar.lambdaJ * vecJ[i]
    end
    for i=(LL-q+1):LL
       grad[i] = 2.0 * plmvar.lambdaH * vecJ[i]
    end 

    vecene = zeros(Float64,q)
    expvecenesunorm = zeros(Float64,q)
    pseudolike = 0.0
    #myrange = Int[tuple((1:site-1)...,(site+1:N)...)...] # concatenate tuples to vector 
    # the cast to Int is relevant: factor 10x in speed
    @inbounds begin 
        for a = 1:M       
            fillvecene!(vecene, vecJ,site,a, q, Z,N)        
            norm = sumexp(vecene)
            expvecenesunorm = exp(vecene .- log(norm))
            #expvecenesunorm = exp(vecene)/norm
            pseudolike -= W[a] * ( vecene[Z[site,a]] - log(norm) )
            offset = 0
            
            
            for i = 1:site-1 #myrange
                @simd for s = 1:q
                    grad[ offset + s + q * ( Z[i,a] - 1 ) ] += W[a] *  expvecenesunorm[s]
                end
                grad[ offset + Z[site,a] + q * ( Z[i,a] - 1 ) ] -= W[a] 
                offset += q2 
            end
	    for i = site+1:N #myrange
                @simd for s = 1:q
                    grad[ offset + s + q * ( Z[i,a] - 1 ) ] += W[a] *  expvecenesunorm[s]
                end
                grad[ offset + Z[site,a] + q * ( Z[i,a] - 1 ) ] -= W[a] 
                offset += q2 
            end

            @simd for s = 1:q 
                grad[ offset + s ] += W[a] *  expvecenesunorm[s] 
            end
	    grad[ offset + Z[site,a] ] -= W[a] 	
        end
    end
    pseudolike += L2norm(vecJ, plmvar)
    return pseudolike 
end

function fillvecene!(vecene::Array{Float64,1}, vecJ::Array{Float64,1},site::Int, a::Int, q::Int, sZ::DenseArray{Int,2},N::Int)
    q2 = q*q
   
    Z = sdata(sZ)
    @inbounds begin
        for l = 1:q
            offset::Int = 0
            scra::Float64 = 0.0
            @simd for i = 1:site-1
                scra += vecJ[offset + l + q * (Z[i,a]-1)] 
            offset += q2 
            end
    	    @simd for i = site+1:N
                scra += vecJ[offset + l + q * (Z[i,a]-1)] 
            offset += q2 
            end
            scra += vecJ[offset + l] # sum H 
            vecene[l] = scra
        end
    end
end

function sumexp(vec::Array{Float64,1})

    mysum = 0.0
    @inbounds begin
        @simd for i=1:length(vec)
            mysum += exp(vec[i])
        end
    end
    return mysum
end

function L2norm(vec::Array{Float64,1}, plmvar::PlmVar)

    LL = length(vec)
    mysum1 = 0.0
    for i=1:(LL-plmvar.q)
        mysum1 += vec[i] * vec[i]
    end
    mysum1 *= plmvar.lambdaJ

    mysum2 = 0.0
    for i=(LL-plmvar.q+1):LL
        mysum2 += vec[i] * vec[i]
    end
    mysum2 *= plmvar.lambdaH
    
    return mysum1+mysum2
end

nothing 
