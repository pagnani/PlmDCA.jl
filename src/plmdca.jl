using GaussDCA
using Optim

immutable type PlmAlg
    epsconv::Float64
    epsgrad::Float64
    epsval::Float64
    maxit::Int
end
immutable type PlmVar
    N::Int
    M::Int
    q::Int    
    q2::Int
    lambdaJ::Float64
    lambdaH::Float64
    Z::Array{Int,2}
    W::Array{Float64,1}
end

function plmdca(filename::String; 
                max_gap_fraction::Real = 0.9, 
                theta = :auto, 
                lambdaJ::Real=0.005, 
                lambdaH::Real=0.01,
                epsconv::Real=0.5,
                epsgrad::Real=1e-5,
                epsval::Real=1e-5,
                maxit::Int=1000,
                method::Symbol=:cg)
    
    W,Z,N,M,q = ReadFasta(filename,max_gap_fraction, theta)
    
    plmalg = PlmAlg(epsconv, epsgrad,epsval,maxit)
    plmvar = PlmVar(N,M,q,q*q,lambdaJ,lambdaH,Z,W)
    
    initvec = zeros(Float64,(plmvar.N - 1) * plmvar.q2 + plmvar.q)
    Jmat = zeros(Float64,(plmvar.N - 1) * plmvar.q2 + plmvar.q, N)
    for site = 1:N
        function f(x::Vector)
            val = PLsiteAndGrad!(x, site, initvec, plmvar)
            return val
        end
        function g!(x::Vector, storage::Vector)
            val = PLsiteAndGrad!(x, site, storage, plmvar)
        end
        @time res = optimize(f, g!, initvec, method=method, show_trace=false, ftol=plmalg.epsval,grtol=plmalg.epsgrad)
        Jmat[:,site] = copy(res.minimum)
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

function GradientDescent(site::Int, plmvar::PlmVar, plmalg::PlmAlg)

    vecJ = zeros(Float64, (plmvar.N - 1) * plmvar.q2 + plmvar.q) 
    grad = zeros(Float64, (plmvar.N - 1) * plmvar.q2 + plmvar.q) 
    it = 1
    normgrad = 10.
    val=0.0
    oldval = val
    deltaval = 1000000.0
    while it <= plmalg.maxit && normgrad > plmalg.epsgrad && deltaval > plmalg.epsval
        val = PLsiteAndGrad!(vecJ, site, grad, plmvar)
        vecJ -= plmalg.epsconv * grad
        normgrad = sum(grad .* grad)
        it += 1
        deltaval = abs(oldval - val)
        oldval = val
    end
    println("it ", it, " |grad| = ", normgrad, " val = ", val, " deltaval = ", deltaval/val)
    #println("it ", it, " |grad| = ", normgrad, " val = ", val, )
    return vecJ
end

function PLsiteAndGrad!(vecJ::Array{Float64,1}, site::Int, grad::Array{Float64,1}, plmvar::PlmVar)

    LL = length(vecJ)
    q2 = plmvar.q2
    q = plmvar.q
    N = plmvar.N
    M = plmvar.M
    Z = plmvar.Z
    W = plmvar.W

    for i=1:LL-q
        grad[i] = 2.0 * plmvar.lambdaJ * vecJ[i]
    end
    for i=(LL-q+1):LL
       grad[i] = 2.0 * plmvar.lambdaH * vecJ[i]
    end 

    vecene = zeros(Float64,q)
    expvecenesunorm = zeros(Float64,q)
    pseudolike = 0.0
    myrange = Int[tuple((1:site-1)...,(site+1:N)...)...] # concatenate tuples to vector 
    # the cast to Int is relevant: factor 10x in speed
    @inbounds begin 
        for a = 1:M       
            fillvecene!(vecene, vecJ, myrange,a, q, plmvar.Z)        
            norm = sumexp(vecene)
            expvecenesunorm = exp(vecene .- log(norm))
            #expvecenesunorm = exp(vecene)/norm
            pseudolike -= W[a] * ( vecene[Z[site,a]] - log(norm) )
            offset = 0
            
            
            for i = myrange
                @simd for s = 1:q
                    grad[ offset + s + q * ( Z[i,a] - 1 ) ] -= W[a] * ( ID(Z[site,a],  s)  -  expvecenesunorm[s])
                end
                offset += q2 
            end
            for s = 1:q 
                grad[ offset + s ] -= W[a] * ( ID(Z[site,a],s) - expvecenesunorm[s] )
            end
        end
    end
    pseudolike += L2norm(vecJ, plmvar)
    return pseudolike 
end

function fillvecene!(vecene::Array{Float64,1}, vecJ::Array{Float64,1},myrange::Array{Int,1}, a::Int, q::Int, Z::Array{Int,2})
    q2 = q*q
    @inbounds begin
        for l = 1:q
            offset = 0
            scra = 0.0
            for i = myrange  # sum J
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
        for i=1:length(vec)
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

function ID(x::Int, y::Int)    
    return (x == y ? 1.0 : 0.0)
end


#a = @parallel vcat for i = 1:N
#           i,f()
#       end
