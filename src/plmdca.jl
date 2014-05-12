using GaussDCA
using Optim
using OptionsMod

immutable type PlmAlg
    method::Symbol
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
                epsconv::Real=0.5,
                epsgrad::Real=1e-5,
                epsval::Real=1e-5,
                maxit::Int=1000,
                method::Symbol=:cg)
    
    W,Z,N,M,q = ReadFasta(filename,max_gap_fraction, theta)

    plmalg = PlmAlg(method,epsconv, epsgrad,epsval,maxit)
    plmvar = PlmVar(N,M,q,q*q,lambdaJ,lambdaH,Z,W)
    if method == :cg2 
        Jmat = MinimizePL2(plmalg, plmvar)
    else
        Jmat = MinimizePL(plmalg, plmvar)
    end

#    return Jmat
#    TransformToMatrix(Jmat,PlmVar)



end
    
# function TransformToMatrix(Jmat::Array{Float64,2}, var::PlmVar)

#     J = zeros(Float64, var.N * var.q, var.N * var.q)
#     h = zeros(Float64, var.N * var.q)



#     for site=1:var.N
#         offseti = 0
#         for i=1:var.N
#             if i != site
#                 for a = 1:var.q
#                     for b = 1:var.q
#                         J[offseti + b]  = vecJ[offsetvecJ + a + var.q*(b-1),site]
#                     end
#                 end
#             end
#             offsetvecJ += var.q2
#             offsetsite += var.q           
#         end
#         offsetsite += 
#     end
# end






end


function MinimizePL2(alg::PlmAlg, var::PlmVar)

    
    x0 = zeros(Float64,(var.N - 1) * var.q2 + var.q)

    Jmat = @parallel hcat for site=1:12#var.N
       function f(g, x::Vector)
            if g === nothing
                g = zeros(Float64, length(x))
            end
            return PLsiteAndGrad!(x, site, g, var)            
        end
        ops = @options display=false fcountmax=alg.maxit tol=alg.epsval
        elapstime = @elapsed  x, fval, fcount, converged = cgdescent(f, x0, ops)
        
        @printf("site = %d\t num-iter = %d\t pl = %.4f\t time = %.4f\n", site, fcount[end],  fval[end], elapstime)
        x
    end 
    return Jmat
end



function MinimizePL(alg::PlmAlg, var::PlmVar)
   
    initvec = zeros(Float64,(var.N - 1) * var.q2 + var.q)

    Jmat = @parallel hcat for site = 1:12#var.N
        function f(x::Vector)
            val = PLsiteAndGrad!(x, site, initvec, var)
            return val
        end
        function g!(x::Vector, storage::Vector)
            val = PLsiteAndGrad!(x, site, storage, var)
        end
        elapstime = @elapsed res = optimize(f, g!, initvec, method=alg.method, show_trace=false, ftol=alg.epsval,grtol=alg.epsgrad)
        @printf("site = %d\t num-iter = %d\t pl = %.4f\t time = %.4f\n", site, res.iterations,  res.f_minimum, elapstime)
        res.minimum
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
    myrange = Int[tuple((1:site-1)...,(site+1:N)...)...] # concatenate tuples to vector 
    # the cast to Int is relevant: factor 10x in speed
    @inbounds begin 
        for a = 1:M       
            fillvecene!(vecene, vecJ, myrange,a, q, Z)        
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

function fillvecene!(vecene::Array{Float64,1}, vecJ::Array{Float64,1},myrange::Array{Int,1}, a::Int, q::Int, sZ::DenseArray{Int,2})
    q2 = q*q
    
    Z = sdata(sZ)
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

nothing 
#a = @parallel vcat for i = 1:N
#           i,f()
#       end
