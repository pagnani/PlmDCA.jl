immutable PlmAlg
    method::Symbol
    verbose::Bool
    epsconv::Float64
    maxit::Int
    boolmask::Union(SharedArray{Bool,2},Nothing)
    function PlmAlg(method,verbose, epsconv, maxit, boolmask)
        if boolmask != nothing 
            sboolmask = SharedArray(Bool,size(boolmask))
            sboolmask[:] = boolmask
            new(method, verbose, epsconv, maxit, sboolmask)
        else
            new(method, verbose, epsconv, maxit, nothing)
        end
    end
end

immutable PlmOut{N}
    pslike::Float64
    Jtensor::Array{Float64,N}
    score::Array{(Int, Int, Float64),1}  
end

immutable PlmVar
    N::Int
    M::Int
    q::Int    
    q2::Int
    gaugecol::Int
    lambdaJ::Float64
    lambdaH::Float64
    Z::SharedArray{Int,2}
    W::SharedArray{Float64,1}
    function PlmVar(N,M,q,q2,gaugecol,lambdaJ, lambdaH, Z,W)
        sZ = SharedArray(Int,size(Z))
        sZ[:] = Z
        sW = SharedArray(Float64,size(W))
        sW[:] = W
        new(N,M,q,q2,gaugecol,lambdaJ, lambdaH, sZ,sW)
    end
end

immutable DecVar{N}
    fracdec::Float64
    fracmax::Float64
    dmask::SharedArray{Bool,N}
    function DecVar(fracdec, fracmax, dmask)
        sdmask = SharedArray(Bool,size(dmask))
        sdmask[:] = dmask 
        new(fracdec, fracmax, sdmask)
    end    
end
