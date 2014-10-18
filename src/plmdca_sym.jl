function plmdcasym(filename::String;
                decimation::Bool=false,
                fracmax::Float64 = 0.3,
                fracdec::Float64 = 0.1,
                remove_dups::Bool = true,
                max_gap_fraction::Real = 0.9, 
                theta = :auto, 
                lambdaJ::Real=0.01, 
                lambdaH::Real=0.01,
                gaugecol::Int=-1,
                epsconv::Real=1.0e-5,
                maxit::Int=1000,
                verbose::Bool=true,
                method::Symbol=:LD_LBFGS)
    
    gaugecol >= 1 && println("Warning: gaugecol not implemented. Proceeding ...")

    W,Z,N,M,q = ReadFasta(filename,max_gap_fraction, theta, remove_dups)

    plmalg = PlmAlg(method,verbose, epsconv ,maxit, nothing)
    plmvar = PlmVar(N,M,q,q*q,gaugecol,lambdaJ,lambdaH,Z,W)

    if decimation == false
        Jmat, pseudolike = MinimizePLSym(plmalg,plmvar) 
    else
        decvar = DecVar{1}(fracdec, fracmax, ones(Bool, binomial(N,2)*q*q+N*q))        
        Jmat, pseudolike = DecimateSym!(plmvar, plmalg, decvar)
    end

    score, Jtens = ComputeScoreSym(Jmat, plmvar)
    return output = PlmOut{3}(pseudolike, Jtens, score)    
end
    
function MinimizePLSym(alg::PlmAlg, var::PlmVar)

    N  = var.N
    q  = var.q
    q2 = var.q2
    
    Nc2 = binomial(N,2)
    LL  = Nc2 * q2  + N * q 

    x0 = zeros(Float64, LL)
    
    function f(x::Vector, g::Vector)
        eltime = @elapsed begin
            g === nothing && (g = zeros(Float64, length(x)))
            pl = PLsiteAndGradSym!(x, g, var)
        end
        alg.verbose && println("time spent in iteration = $eltime . PL = $pl")
        return pl
    end

    opt = Opt(alg.method, length(x0))
    
    ftol_abs!(opt, alg.epsconv)
    maxeval!(opt, alg.maxit)
    min_objective!(opt, f)
    elapstime = @elapsed  (minf, minx, ret) = optimize(opt, x0)
    @printf("pl = %.4f\t time = %.4f\t exit status = ", minf, elapstime)
    println(ret)
    
    return minx, minf
end

function ComputeScoreSym(Jvec::Array{Float64,1}, var::PlmVar)


    LL = length(Jvec)
    N=var.N
    q=var.q
    Nc2 = binomial(N,2)

    Jtens=reshape(Jvec[1:LL-N*q],q,q,Nc2)
    J = zeros(q,q,Nc2)

    for l=1:Nc2
        J[:,:,l] = Jtens[:,:,l] - repmat(mean(Jtens[:,:,l],1),q,1)-repmat(mean(Jtens[:,:,l],2),1,q) .+ mean(Jtens[:,:,l])
    end

    FN = zeros(Float64, N,N)
    l = 1
    for i=1:N-1
        for j=i+1:N
            FN[i,j] = vecnorm(J[:,:,l],2)
            FN[j,i] = FN[i,j]
            l+=1
        end
    end
    FN = GaussDCA.correct_APC(FN)  
    score = GaussDCA.compute_ranking(FN)
    return score, Jtens
end

function PLsiteAndGradSym!(vecJ::Array{Float64,1}, grad::Array{Float64,1}, plmvar::PlmVar)

    LL = length(vecJ)
    q2 = plmvar.q2
    q = plmvar.q
    N = plmvar.N
    M = plmvar.M
    Z = plmvar.Z
    W = plmvar.W
    lambdaJ = plmvar.lambdaJ
    lambdaH = plmvar.lambdaH

    for i=1:LL-N*q
        grad[i] = 2.0 * vecJ[i] * lambdaJ
    end
    for i=(LL-N*q + 1):LL
        grad[i] = 2.0 * vecJ[i] * lambdaH
    end
    pseudolike = 0
    for a = 1:M         
        pseudolike += ComputePatternPLSym!(grad, vecJ, Z[:,a], W[a], N, q, q2)
    end
    
    pseudolike += L2norm_sym(vecJ, plmvar)
    return pseudolike 
end

function ComputePatternPLSym!(grad::Array{Float64,1}, vecJ::Array{Float64,1}, Z::Array{Int,1}, Wa::Float64, N::Int, q::Int, q2::Int)
    
    vecene = zeros(Float64,q)
    expvecenesunorm = zeros(Float64,q)
    pseudolike = 0
    offset = mygetindex(N-1, N, q, q, N, q, q2) 
    @inbounds for site=1:N    # site < i 
        fillvecenesym!(vecene, vecJ, Z, site, q,N)        
        norm = sumexp(vecene)
        expvecenesunorm = exp(vecene .- log(norm))
        pseudolike -= Wa * ( vecene[Z[site]] - log(norm) )
	for i = 1:(site-1)
            for s = 1:q
                grad[ mygetindex(i, site, Z[i], s, N, q, q2) ] += 0.5 * Wa * expvecenesunorm[s]
            end
            grad[ mygetindex(i, site , Z[i], Z[site],  N,q,q2)] -= 0.5 * Wa
        end
	for i = (site+1):N 
            for s = 1:q
                grad[ mygetindex(site, i , s,  Z[i], N,q,q2) ] += 0.5 * Wa * expvecenesunorm[s]
            end
            grad[ mygetindex(site, i , Z[site], Z[i], N,q,q2)] -= 0.5* Wa
        end
        @simd for s = 1:q 
            grad[ offset + s ] += Wa *  expvecenesunorm[s] 
        end
	grad[ offset + Z[site] ] -= Wa
 	offset += q
    end
    return pseudolike
end


function fillvecenesym!(vecene::Array{Float64,1}, vecJ::Array{Float64,1}, Z::Array{Int64,1}, site::Int, q::Int ,N::Int)
    q2 = q*q
#    Z = sdata(sZ)

    @inbounds begin
        for l = 1:q
            offset::Int = 0
            scra::Float64 = 0.0

            for i=1:1:site-1
                scra += vecJ[ mygetindex(i, site, Z[i], l,  N, q, q2)]
            end
    	    for i = site+1:N
                scra += vecJ[ mygetindex(site, i, l, Z[i], N, q, q2)]
            end # End sum_i \neq site J
            scra *= 0.5 
            offset = mygetindex(N-1, N, q, q, N, q, q2)  + ( site - 1) * q  # last J element + (site-1)*q
#            println(length(vecJ), " offset ", offset) 
            scra += vecJ[offset + l] # sum H 
            vecene[l] = scra
        end
    end
end

function L2norm_sym(vec::Array{Float64,1}, var::PlmVar)

    q = var.q    
    N = var.N
    lambdaJ = var.lambdaJ
    lambdaH = var.lambdaH

    LL = length(vec)


    mysum1 = 0.0
    @inbounds @simd for i=1:(LL-N*q)
        mysum1 += vec[i] * vec[i]
    end
    mysum1 *= lambdaJ

    mysum2 = 0.0
    @inbounds @simd for i=(LL - N*q + 1):LL
        mysum2 += vec[i] * vec[i]
    end
    mysum2 *= lambdaH
    
    return mysum1+mysum2
end


function mygetindex( i::Int, j::Int, coli::Int, colj::Int, N::Int, q::Int, q2::Int)        
    offset_i = ( (i-1) * N  - ( (i * ( i -1 ) ) >> 1 ) ) * q2 # (i-1) N q2 + i (i-1) q2 / 2  
    offset_j = (j - i - 1 ) * q2
    return offset_i + offset_j + coli + q * (colj - 1)
end


nothing
