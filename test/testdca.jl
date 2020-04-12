module TestDCA
using Test, PottsGauge

function energy(x,J,h)
    q,q,N,N = size(J)
    eJ = 0.0
    @inbounds @simd for i in 1:N-1
        for j in i+1:N
            eJ -= J[x[i],x[j],i,j]
        end
    end
    eh = 0.0
    @inbounds @simd for i in 1:N-1
        eh -= h[x[i],i]
    end
    return eJ + eh
end

function complete_dataset(N,q)
    Z = zeros(Int,N,q^N)
    bounds = ntuple(x->q,N)
    ctr =0;
    @inbounds for i in CartesianIndices(bounds)
        ctr += 1
        for j in 1:N
            Z[j,ctr] = i[j]
        end
    end
    Z
end

function generateWZJh(N,q)
    Jasym = randn(q,q,N,N)
    h = randn(q,N)
    J = 0.5*(permutedims(Jasym,[2,1,4,3])+Jasym)
    for i in 1:N
        for a in 1:q
            for b in 1:q
                J[a,b,i,i] = 0.0
            end
        end
    end
    Z = complete_dataset(N,q)
    W = [exp(-energy(Z[:,i],J,h)) for i in 1:q^N]
    W.= W ./sum(W)
    return W,Z,J,h
end

function myplmDCA(Z::Array{T,2}, W::Vector{S};
    method::Symbol=:LD_LBFGS,
    verbose::Bool = true,
    min_separation::Integer=1,
    epsconv::Real = 1e-5,
    maxit::Integer = 1000,
    boolmask = nothing,
    gaugecol::Int = -1,
    lambdaJ::Real = 0.01,
    lambdaH::Real = 0.01) where T<:Integer where S<:AbstractFloat

    N, M = size(Z)
    q = Int(maximum(Z))
    q > 32 && error("parameter q=$q is too big (max 32 is allowed)")
    M == length(W) || throw(DimensionMismatch("incompatible length W"))
    plmalg = PlmDCA.PlmAlg(method, verbose, epsconv , maxit, boolmask)
    normedW = W./sum(W) # plmDCA expects normalized W !!!!
    plmvar = PlmDCA.PlmVar(N,M,q,q*q,gaugecol,lambdaJ,lambdaH,round.(Int,Z),normedW)
    Jmat, pslike = PlmDCA.MinimizePLAsym(plmalg,plmvar)
    plmscore, FNAPC, Jtensor, htensor = PlmDCA.ComputeScore(Jmat, plmvar, min_separation)
    return plmscore,Jtensor,htensor
end

function testDCA(N,q)
    W,Z,J,h = generateWZJh(N,q)
    plmscore,Jplm,hplm = myplmDCA(Z,W,lambdaJ=0.0,lambdaH=0.0,epsconv=1e-20,verbose=false)
    Jz,hz = PottsGauge.gauge(J,h,ZeroSumGauge())
    Jplmz,hplmz = PottsGauge.gauge(Jplm,hplm,ZeroSumGauge())
    return sum(abs2,Jz-Jplmz)<1e-6
end

@test testDCA(4,2)
@test testDCA(6,2)
@test testDCA(4,3)
printstyled("All TestDCA passed!\n",color=:green,bold=true)
end
