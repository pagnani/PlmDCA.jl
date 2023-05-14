function plmdca_asym(Z::Array{T,2},W::Vector{Float64};
    remove_dups::Bool=true,
    min_separation::Int=1,
    theta=:auto,
    lambdaJ::Real=0.01,
    lambdaH::Real=0.01,
    epsconv::Real=1.0e-5,
    maxit::Int=1000,
    verbose::Bool=true,
    method::Symbol=:LD_LBFGS
) where {T <: Integer}

    all(x -> x > 0, W) || throw(DomainError("vector W should normalized and with all positive elements"))
    isapprox(sum(W), 1) || throw(DomainError("sum(W) ≠ 1. Consider normalizing the vector W"))
    N, M = size(Z)
    M = length(W)
    q = Int(maximum(Z))
    plmalg = PlmAlg(method, verbose, epsconv, maxit)
    plmvar = PlmVar(N, M, q, lambdaJ, lambdaH, Z, W)
    Jmat, pslike = minimize_pl_asym(plmalg, plmvar)
    score, _, Jtensor, htensor =  compute_score(Jmat, plmvar, min_separation)
    return PlmOut(pslike, Jtensor, htensor, score)
end

plmdca(Z,W;kwds...) = plmdca_asym(Z, W;kwds...)

function plmdca_asym(filename::String;
    theta::Union{Symbol,Real}=:auto,
    max_gap_fraction::Real=0.9,
    remove_dups::Bool=true,
    kwds...
)
    time = @elapsed W, Z, N, M, q = read_fasta(filename, max_gap_fraction, theta, remove_dups)
    println("preprocessing took $time seconds")
    plmdca_asym(Z, W; kwds...)
end

plmdca(filename::String; kwds...) = plmdca_asym(filename; kwds...)

function minimize_pl_asym(alg::PlmAlg, var::PlmVar)

    LL = (var.N - 1) * var.q2 + var.q
    x0 = zeros(Float64, LL)
    g = zeros(Float64, LL)
    vecps = Vector{Float64}(undef, var.N)
    Jmat = zeros(LL,var.N)
    Threads.@threads for site = 1:var.N
        opt = Opt(alg.method, length(x0))
        ftol_abs!(opt, alg.epsconv)
        xtol_rel!(opt, alg.epsconv)
        xtol_abs!(opt, alg.epsconv)
        ftol_rel!(opt, alg.epsconv)
        maxeval!(opt, alg.maxit)
        min_objective!(opt, (x, g) -> optimfunwrapper(x, g, site, var))
        elapstime = @elapsed (minf, minx, ret) = optimize(opt, x0)
        alg.verbose && @printf("site = %d\t pl = %.4f\t time = %.4f\t", site, minf, elapstime)
        alg.verbose && println("exit status = $ret")
        vecps[site] = minf
        Jmat[:,site] .= minx
    end
    return Jmat, vecps
end

function pl_site_grad!(x::Vector{Float64}, grad::Vector{Float64}, site::Int, plmvar::PlmVar)
    LL = length(x)
    q2 = plmvar.q2
    q = plmvar.q
    N = plmvar.N
    M = plmvar.M
    Z = plmvar.Z
    W = plmvar.W
    IdxZ = plmvar.IdxZ
    for i = 1:LL - q
        grad[i] = 2.0 * plmvar.lambdaJ  * x[i]
    end
    for i = (LL - q + 1):LL
        grad[i] = 4.0 * plmvar.lambdaH * x[i]
    end
    pseudolike = 0.0
    vecene = zeros(Float64, q)
    lnorm = 0.0
    expvecenesumnorm = zeros(Float64, q)
    @inbounds for m in 1:M
        izm = view(IdxZ, :, m)
        zsm = Z[site,m]
        fillvecene!(vecene, x, site, izm, q, N)
        lnorm = logsumexp(vecene)
        expvecenesumnorm .= @. exp(vecene - lnorm)
        pseudolike -= W[m] * (vecene[ zsm ] - lnorm)
        @turbo for i in 1:site - 1
        for s = 1:q
                grad[ izm[i] + s ] += W[m] * expvecenesumnorm[s]
            end
            grad[ izm[i] + zsm ] -= W[m]
        end
        @turbo for i = site + 1:N
            for s = 1:q
                grad[ izm[i] - q2 + s ] += W[m] *  expvecenesumnorm[s]
            end
            grad[ izm[i] - q2 + zsm ] -= W[m]
        end
        @turbo for s = 1:q
            grad[ (N - 1) * q2 + s ] += W[m] * expvecenesumnorm[s]
        end
        grad[ (N - 1) * q2 + zsm ] -= W[m]
    end
    pseudolike += l2norm_asym(x, plmvar)
end


# Energy filling

function fillvecene!(vecene::Vector{Float64}, x::Vector{Float64}, site::Int, IdxSeq::AbstractArray{Int,1}, q::Int, N::Int)
	q2 = q * q
    @inbounds for l in 1:q
        scra::Float64 = 0.0
        @turbo for i = 1:site - 1 # Begin sum_i \neq site J
            scra += x[IdxSeq[i] + l]
        end
        # skipping sum over residue site
        @turbo for i = site + 1:N
            scra += x[IdxSeq[i] - q2 + l]
        end # End sum_i \neq site J
        scra +=  x[(N - 1) * q2 + l] # sum H
        vecene[l] = scra
    end
end


function logsumexp(X::Vector{Float64})
    u = maximum(X)
    isfinite(u) || return float(u)
    return u + log(sum(x -> exp(x - u), X))
end


function l2norm_asym(vec::Array{Float64,1}, plmvar::PlmVar)
    q = plmvar.q
    N = plmvar.N
    lambdaJ = plmvar.lambdaJ
    lambdaH = plmvar.lambdaH

    LL = length(vec)
    mysum1 = 0.0
    @inbounds @turbo for i in 1:(LL - q)
        mysum1 += vec[i] * vec[i]
    end
    mysum1 *= lambdaJ
    mysum2 = 0.0
    @inbounds @turbo for i in (LL - q + 1):LL
        mysum2 += vec[i] * vec[i]
    end
    mysum2 *= 2lambdaH
    return mysum1 + mysum2
end



function compute_score(Jmat::Array{Float64,2}, var::PlmVar, min_separation::Int)
    q = var.q
    N = var.N
    JJ = reshape(Jmat[1:end - q,:], q, q, N - 1, N)
    Jtemp1 = zeros(q, q, Int(N * (N - 1) / 2))
    Jtemp2 = zeros(q, q, Int(N * (N - 1) / 2))
    l = 1
    for i = 1:(N - 1)
        for j = (i + 1):N
            Jtemp1[:,:,l] = JJ[:,:,j - 1,i] # J_ij as estimated from from g_i.
            Jtemp2[:,:,l] = JJ[:,:,i,j]' # J_ij as estimated from from g_j.
            l = l + 1
        end
    end

    hplm = fill(0.0, q, N)
    for i in 1:N
        hplm[:,i] .= Jmat[end - q + 1:end,i]
    end
    Jtensor1 = inflate_matrix(Jtemp1, N)
    Jtensor2 = inflate_matrix(Jtemp2, N)
    Jplm = (Jtensor1 + Jtensor2) / 2 # for the energy I do not want to gauge
    ctr = 0
    for i in 1:N - 1
        for j in i + 1:N
            ctr += 1
            Jtensor1[:,:,i,j] = Jtemp1[:,:,ctr] - repeat(mean(Jtemp1[:,:,ctr], dims=1), q, 1) - repeat(mean(Jtemp1[:,:,ctr], dims=2), 1, q) .+ mean(Jtemp1[:,:,ctr])
            Jtensor1[:,:,j,i] = Jtensor1[:,:,i,j]'
            Jtensor2[:,:,i,j] = Jtemp2[:,:,ctr] - repeat(mean(Jtemp2[:,:,ctr], dims=1), q, 1) - repeat(mean(Jtemp2[:,:,ctr], dims=2), 1, q) .+ mean(Jtemp2[:,:,ctr])
            Jtensor2[:,:,j,i] = Jtensor2[:,:,i,j]'
        end
    end # zerosumgauge the different tensors
    Jtensor = (Jtensor1 + Jtensor2) / 2
    FN = compute_APC(Jtensor, N, q)
    score = compute_ranking(FN, min_separation)
    return score, FN, Jplm, hplm
end
