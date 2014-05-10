using GaussDCA

function plmdca(nomefile::String; 
                max_gap_fraction::Real = 0.9 , 
                theta = :auto, 
                lambdaJ::Real=0.005, 
                lambdaH::Real=0.01,
                epsconv::Real=0.1,
                maxit::Int=10,
                epsgrad::Real=1e-5,
                epsval::Real=1e-5)
    
    Z = GaussDCA.read_fasta_alignment(nomefile, max_gap_fraction)
    N, M = size(Z)
    q = int(maximum(Z))

    q > 32 && error("parameter q=$q is too big (max 31 is allowed)")
    _, _, Meff, W = GaussDCA.compute_new_frequencies(Z, theta)
    W  ./= Meff  
    Zint = int( Z )
    

    
    vecJ = zeros(Float64, (N-1) * q*q + q) 
    grad = zeros(Float64, (N-1) * q*q + q) 
    it = 1
    
    # function PL(vecJ)
    #     PLsite(vecJ, site, Zint, W, q, lambdaH, lambdaJ )
    # end
    # function gradPL!(vecJ)
    #     gradPLsite!(vecJ, site, Zint, W, q, lambdaH, lambdaJ, grad)
    # end
    normgrad = 10
    while true && it <= maxit && normgrad > epsgrad
        site = 1
        val = PLsiteAndGrad!(vecJ, site, Zint, W, q, lambdaH, lambdaJ, grad)
        vecJ -= epsconv * grad
        normgrad = sum(grad .* grad)
        println("it ", it, " |grad| = ", normgrad, " val = ", val)
        it += 1
    end
    return Z, vecJ
end

# function converge()

#     normgrad = 10
#     while true && it <= maxit && normgrad > epsgrad
#         site = 1
#         val = PLsiteAndGrad!(vecJ, site, Zint, W, q, lambdaH, lambdaJ, grad)
#         vecJ -= epsconv * grad
#         normgrad = sum(grad .* grad)
#         println("it ", it, " |grad| = ", normgrad, " val = ", val)
#         it += 1
#     end
# end

function PLsiteAndGrad!(vecJ::Array{Float64,1}, site::Int, Z::Array{Int,2}, W::Array{Float64,1}, q::Int, lambdaH::Float64, lambdaJ::Float64, grad::Array{Float64,1} )

    LL = length(vecJ)
    M = length(W)
    N = size(Z,1)

    q2 = q*q
    
    for i=1:LL-q
        grad[i] = 2.0 * lambdaJ * vecJ[i]
    end
    for i=(LL-q+1):LL
       grad[i] = 2.0 * lambdaH * vecJ[i]
    end 

    vecene = zeros(Float64,q)
    expvecenesunorm = zeros(Float64,q)
    pseudolike = 0.0
    myrange = Int[tuple((1:site-1)...,(site+1:N)...)...] # concatenate tuples to vector 
    # the cast to Int is relevant: factor 10x in speed
    @inbounds begin 
        for a = 1:M       
            fillvecene!(vecJ, Z, a, myrange, N, q, vecene)        
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
#    println(pseudolike)
    pseudolike += L2norm(vecJ, lambdaJ, lambdaH, q)
    return pseudolike 
end

function fillvecene!(vecJ::Array{Float64,1}, Z::Array{Int,2}, a::Int, myrange::Array{Int,1}, N::Int, q::Int, vecene::Array{Float64,1})

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

function L2norm(vec::Array{Float64,1}, lambdaJ::Float64, lambdaH::Float64,q::Int)

    LL = length(vec)
    mysum1 = 0.0
    for i=1:(LL-q)
        mysum1 += vec[i] * vec[i]
    end
    mysum1 *= lambdaJ

    mysum2 = 0.0
    for i=(LL-q+1):LL
        mysum2 += vec[i] * vec[i]
    end
    mysum2 *= lambdaH

    return mysum1+mysum2
end

function ID(x::Int, y::Int)    
    return (x == y ? 1.0 : 0.0)
end


#a = @parallel vcat for i = 1:N
#           i,f()
#       end
