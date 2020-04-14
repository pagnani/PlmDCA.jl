module TestDCA
using Test, PlmDCA, PottsGauge

function energy(x,J,h)
    q,q,N,N = size(J)
    eJ = 0.0
    @inbounds @simd for i in 1:N-1
        for j in i+1:N
            eJ -= J[x[i],x[j],i,j]
        end
    end
    eh = 0.0
    @inbounds @simd for i in 1:N
        eh -= h[x[i],i]
    end
    return eJ + eh
end

function complete_dataset(N,q)
    Z = zeros(Int,N,q^N)
    bounds = ntuple(x->q,N)
    ctr = 0
    @inbounds for i in CartesianIndices(bounds)
        ctr += 1
        for j in 1:N
            Z[j,ctr] = i[j]
        end
    end
    Z
end

function generateWZJh(N,q)
    Jasym = rand(q,q,N,N)
    h = rand(q,N)
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

function testDCA(N,q;
                 verbose::Bool=false,
                 epsconv::Real=1e20, gauge = ZeroSumGauge(),
                 lambdaJ::Real = 0.0,
                 lambdaH=0.0,
                 asym::Bool=true,
                 maxit::Integer=1000,
                 method::Symbol=:LD_LBFGS)
   
    W,Z,J,h = generateWZJh(N,q)

    func_dca = asym ? plmdca_asym : plmdca_sym
    
    resplm = func_dca(Z,W,lambdaJ=lambdaJ,lambdaH=lambdaH,epsconv=epsconv,verbose=verbose,maxit=maxit,method=method)
   
    Jplm,hplm = resplm.Jtensor,resplm.htensor
    Jz,hz = PottsGauge.gauge(J,h,gauge)
    Jplmz,hplmz = PottsGauge.gauge(Jplm,hplm,gauge)
#    return Jplmz,hplmz,Jz,hz,Jplm,hplm,J,h
    return sum(abs2,Jz-Jplmz)<1e-6
end
for tf in (true, false) # asymmetric (true) and symmetric (false
    @test testDCA(4,2,lambdaJ=1e-7,epsconv=1e-20,asym=tf,verbose=true)
    @test testDCA(6,2,lambdaJ=1e-7,epsconv=1e-20,asym=tf,verbose=true)
    @test testDCA(4,3,lambdaJ=1e-7,epsconv=1e-20,asym=tf,verbose=true)
end
printstyled("All TestDCA passed!\n",color=:light_green,bold=true)
end
