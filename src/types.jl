struct PlmAlg
    method::Symbol
    verbose::Bool
    epsconv::Float64
    maxit::Int
end

struct PlmOut
    pslike::Union{Vector{Float64},Float64}
    Jtensor::Array{Float64,4}
    htensor::Array{Float64,2}
    score::Array{Tuple{Int, Int, Float64},1}  
end

struct PlmVar
    N::Int
    M::Int
    q::Int    
    q2::Int
    lambdaJ::Float64
    lambdaH::Float64
    Z::SharedArray{Int,2}
    W::SharedArray{Float64,1}
    function PlmVar(N,M,q,q2,lambdaJ, lambdaH, Z,W)
        sZ = SharedArray{Int}(size(Z))
        sZ[:] = Z
        sW = SharedArray{Float64}(size(W))
        sW[:] = W
        new(N,M,q,q2,lambdaJ, lambdaH, sZ,sW)
    end
end

struct DecVar{N} 
    fracdec::Float64
    fracmax::Float64
    blockdecimate::Bool
    dmask::SharedArray{Bool,N}
    function DecVar{N}(fracdec, fracmax, blockdecimate, dmask) where N
        sdmask = SharedArray{Bool}(size(dmask))
        sdmask[:] = dmask 
        new(fracdec, fracmax, blockdecimate, sdmask)
    end
end
