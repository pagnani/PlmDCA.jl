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
    Z::Array{Int,2}
    W::Array{Float64,1}
    IdxZ::Array{Int,2} #partial index computation for speed up energy calculation

    function PlmVar(N,M,q,lambdaJ, lambdaH, Z,W)
        IdxZ = Array{Int,2}(undef, N, M)
        q2=q*q
        for i in 1:M
            for j in 1:N
                IdxZ[j,i] = (j-1) * q2 + q * (Z[j,i] - 1)
            end
        end
        new(N,M,q,q2,lambdaJ, lambdaH, Z, W, IdxZ)
    end
end


# struct DecVar{N} 
#     fracdec::Float64
#     fracmax::Float64
#     blockdecimate::Bool
#     dmask::SharedArray{Bool,N}
#     function DecVar{N}(fracdec, fracmax, blockdecimate, dmask) where N
#         sdmask = SharedArray{Bool}(size(dmask))
#         sdmask[:] = dmask 
#         new(fracdec, fracmax, blockdecimate, sdmask)
#     end
# end
