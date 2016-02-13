function ComputeScore(Jmat::Array{Float64,2}, var::PlmVar, min_separation::Int)

    q = var.q
    N = var.N

    JJ=reshape(Jmat[1:end-q,:], q,q,N-1,N)
@compat    Jtemp1=zeros( q,q,Int(N*(N-1)/2))
@compat    Jtemp2=zeros( q,q,Int(N*(N-1)/2))
    l = 1

    for i=1:(N-1)
        for j=(i+1):N
            Jtemp1[:,:,l]=JJ[:,:,j-1,i]; #J_ij as estimated from from g_i.
            Jtemp2[:,:,l]=JJ[:,:,i,j]'; #J_ij as estimated from from g_j.
            l=l+1;
        end
    end
    
    Jtensor = zeros(q,q,N,N)
    l = 1
    for i = 1:N-1
        for j=i+1:N
            Jtensor[:,:,i,j] = Jtemp1[:,:,l]
            Jtensor[:,:,j,i] = Jtemp2[:,:,l]'
            l += 1
        end
    end

    ASFN = zeros(N,N)
    for i=1:N,j=1:N 
        i!=j && (ASFN[i,j] =sum(Jtensor[:,:,i,j].^2)) 
    end

@compat    J1=zeros(q,q,Int(N*(N-1)/2))
@compat    J2=zeros(q,q,Int(N*(N-1)/2))

@compat    for l=1:Int(N*(N-1)/2)
        J1[:,:,l] = Jtemp1[:,:,l]-repmat(mean(Jtemp1[:,:,l],1),q,1)-repmat(mean(Jtemp1[:,:,l],2),1,q) .+ mean(Jtemp1[:,:,l])
        J2[:,:,l] = Jtemp2[:,:,l]-repmat(mean(Jtemp2[:,:,l],1),q,1)-repmat(mean(Jtemp2[:,:,l],2),1,q) .+ mean(Jtemp2[:,:,l])
    end
    J = 0.5 * ( J1 + J2 )



    FN = zeros(Float64, N,N)
    l = 1

    for i=1:N-1
        for j=i+1:N
            FN[i,j] = vecnorm(J[1:q-1,1:q-1,l],2)
#            FN[i,j] = vecnorm(J[:,:,l],2)
            FN[j,i] =FN[i,j]
            l+=1
        end
    end
    FN=GaussDCA.correct_APC(FN)
    score = GaussDCA.compute_ranking(FN,min_separation)
    return score, FN, Jtensor
end


function ReadFasta(filename::AbstractString,max_gap_fraction::Real, theta::Any, remove_dups::Bool)

    Z = GaussDCA.read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups
        Z, _ = GaussDCA.remove_duplicate_seqs(Z)
    end


    N, M = size(Z)
@compat    q = round(Int,maximum(Z))
    
    q > 32 && error("parameter q=$q is too big (max 31 is allowed)")
    W , Meff = GaussDCA.compute_weights(Z,q,theta)
    scale!(W, 1.0/Meff)
@compat    Zint=round(Int,Z)
    return W, Zint,N,M,q
end

function sumexp(vec::Array{Float64,1})    
    mysum = 0.0
    @inbounds @simd for i=1:length(vec)
        mysum += exp(vec[i])
    end
    return mysum
end



function Base.deepcopy_internal(x::DecVar, d::ObjectIdDict)
    haskey(d, x) && return d[x]
    dmask_c = Base.deepcopy_internal(x.dmask, d)
    xc = DecVar(x.fracdec, x.fracmax, dmask_c)
    d[x] = xc
    return xc
end
