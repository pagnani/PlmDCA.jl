function optimfunwrapper(x::Vector, g::Vector, site, var)
    g === nothing && (g = zeros(Float64, length(x)))
    return PLsiteAndGrad!(x, g, site,  var)
end

function optimfunwrapper(x::Vector, g::Vector, var)
    g === nothing && (g = zeros(Float64, length(x)))
    return PLsiteAndGradSym!(x, g, var)
end


function ComputeScore(Jmat::Array{Float64,2}, var::PlmVar, min_separation::Int)

    q = var.q
    N = var.N

    JJ=reshape(Jmat[1:end-q,:], q,q,N-1,N)
    Jtemp1=zeros( q,q,Int(N*(N-1)/2))
    Jtemp2=zeros( q,q,Int(N*(N-1)/2))
    l = 1
    for i=1:(N-1)
        for j=(i+1):N
            Jtemp1[:,:,l]=JJ[:,:,j-1,i]; #J_ij as estimated from from g_i.
            Jtemp2[:,:,l]=JJ[:,:,i,j]'; #J_ij as estimated from from g_j.
            l=l+1;
        end
    end




    hplm = fill(0.0, q,N)
    for i in 1:N
        hplm[:,i] .= Jmat[end-q+1:end,i]
    end

    Jtensor1 = inflate_matrix(Jtemp1,N)
    Jtensor2 = inflate_matrix(Jtemp2,N)
    Jplm = (Jtensor1 + Jtensor2)/2 # for the energy I do not want to gauge

    ctr = 0
    for i in 1:N-1
        for j in i+1:N
            ctr += 1
            Jtensor1[:,:,i,j] = Jtemp1[:,:,ctr]-repeat(mean(Jtemp1[:,:,ctr],dims=1),q,1)-repeat(mean(Jtemp1[:,:,ctr],dims=2),1,q) .+ mean(Jtemp1[:,:,ctr])
            Jtensor1[:,:,j,i] = Jtensor1[:,:,i,j]'
            Jtensor2[:,:,i,j] = Jtemp2[:,:,ctr]-repeat(mean(Jtemp2[:,:,ctr],dims=1),q,1)-repeat(mean(Jtemp2[:,:,ctr],dims=2),1,q) .+ mean(Jtemp2[:,:,ctr])
            Jtensor2[:,:,j,i] = Jtensor2[:,:,i,j]'
        end
    end # zerosumgauge the different tensors

    Jtensor = (Jtensor1 + Jtensor2)/2

    FN = compute_APC(Jtensor,N,q)
    score = GaussDCA.compute_ranking(FN,min_separation)

    return score, FN, Jplm, hplm
end


function compute_APC(J::Array{Float64,4},N,q)
    FN = fill(0.0, N,N)
    for i=1:N-1
        for j=i+1:N
            FN[i,j] = norm(J[1:q-1,1:q-1,i,j],2)
            FN[j,i] =FN[i,j]
        end
    end
    FN=GaussDCA.correct_APC(FN)
    return FN
end


function inflate_matrix(J::Array{Float64,3},N)
    q,q,NN = size(J)

    @assert (N*(N-1))>>1 == NN

    Jt = zeros(q,q,N,N)
    ctr = 0
    for i in 1:N-1
        for j in i+1:N
            ctr += 1
            Jt[:,:,i,j] = J[:,:,ctr]
            Jt[:,:,j,i] = J[:,:,ctr]'
        end
    end
    return Jt
end




function ReadFasta(filename::AbstractString,max_gap_fraction::Real, theta::Any, remove_dups::Bool)

    Z = GaussDCA.read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups
        Z, _ = GaussDCA.remove_duplicate_seqs(Z)
    end


    N, M = size(Z)
    q = round(Int,maximum(Z))

    q > 32 && error("parameter q=$q is too big (max 31 is allowed)")
    W , Meff = GaussDCA.compute_weights(Z,q,theta)
    rmul!(W, 1.0/Meff)
    Zint=round.(Int,Z)
    return W, Zint,N,M,q
end

function sumexp(vec::Array{Float64,1})
    mysum = 0.0
    @inbounds @simd for i=1:length(vec)
        mysum += exp(vec[i])
    end
    return mysum
end
