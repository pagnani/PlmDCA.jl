function DecimateSym!(var::PlmVar, alg::PlmAlg, dec::DecVar{1})

    Jvec,pseudolike = MinimizePLSym(alg, var)
    Nc2 = binomial(var.N,2)
    blockdecimate = dec.blockdecimate
    LL = blockdecimate ? Nc2 : Nc2 * var.q2
    iteration = 1
    while true        
        SetAllCouplingToZero!(Jvec,  var, dec)
        nzeros = sum(!,dec.dmask) # number of couplings to be set to zero
        fraczeros = nzeros / LL 
        println("# set to zero ", nzeros ," over ", LL ," equal to ", 100.0*fraczeros," % of the couplings. Converging ...")
        pseudolike, Jvec = MinimizePLSym(alg, var, dec, Jvec, x0=Jvec)
        println("# Converged!")
#        push!(results,(iteration,1.0-fraczeros, fmin,fminloc))
        fraczeros >= 1.0 && break
        fraczeros >= dec.fracmax && break
        iteration += 1
    end
    return Jvec, pseudolike
end

function SetAllCouplingToZero!(vecJ::Array{Float64,1}, var::PlmVar, dec::DecVar{1})

    N = var.N
    q  = var.q
    q2 = var.q2
    blockdecimate = dec.blockdecimate

    if blockdecimate 
        LL = binomial(N,2) 
    else
        LL = length(vecJ)
        LL == (binomial(N,2)*q2 + N*q) || error("LL = $LL should be ", binomial(N,2)*q2 + N*q)
    end

    
    vecval = Tuple{Int,Float64}[]
    sizehint!(vecval,LL)
    if blockdecimate   # compute the FN of the block 
        for i=1:LL
            if dec.dmask[i]
                fnvecJ = 0.0
                for j in ((i-1)*q2 + 1):i*q2
                    fnvecJ += vecJ[j]^2
                end
                push!(vecval,(i,fnvecJ))
            end
        end
    else
        for i=1:LL-N*q
            dec.dmask[i] && push!(vecval,(i, vecJ[i]^2))
        end
    end
    
    lvec = length(vecval)
    idxperm = sortperm(Float64[ vecval[i][2] for i=1:lvec]) 
    
    ndec = blockdecimate ? round(Int, dec.fracdec * (LL-N*q)) : round(Int, dec.fracdec*LL)
    
    for i=1:ndec
        site = vecval[idxperm[i]][1]
        dec.dmask[site] = false
    end    
end


function MinimizePLSym(alg::PlmAlg, var::PlmVar, dec::DecVar{1}, Jvec::Array{Float64,1}; x0 = zeros(Float64,length(Jvec)))
    
    N  = var.N
    q  = var.q
    q2 = var.q2
    Nc2 = binomial(N,2)
    
    opt = Opt(alg.method, length(x0))
    ftol_abs!(opt, alg.epsconv)
    lb,ub = ComputeConstraintsSym!(dec, var.q, var.N, x0)
    lower_bounds!(opt, lb)
    upper_bounds!(opt, ub)

    maxeval!(opt, alg.maxit)
    min_objective!(opt, (x,g)->optimfunwrapper(x,g,var))
    elapstime = @elapsed  (minf, minx, ret) = optimize(opt, x0)

    alg.verbose && @printf("pl = %.4f\t time = %.4f\n",  minf, elapstime)
    return minf,minx
end 

function ComputeConstraintsSym!(decvar::DecVar{1}, q::Int, N::Int, x0::Array{Float64,1} )
    
    q2 = q*q    
    dmask = decvar.dmask
    blockdecimate = decvar.blockdecimate
    
    Nc2 = binomial(N,2)
    lb = -Inf * ones(Float64, Nc2*q2 + N*q)
    ub =  Inf * ones(Float64, Nc2*q2 + N*q)
    LL = length(x0)

    if blockdecimate
        for i=1:Nc2
            if dmask[i] == false
                for j in ((i-1)*q2 + 1):i*q2
                    lb[ j ] = -1e-12
                    ub[ j ] = 1e-12
                    x0[ j ] = 0.0
                end
            else
                for j in ((i-1)*q2 + 1):i*q2
                    x0[ j ] += (1.0 - 2.0*rand())*1e-3
                end
            end
        end
    else
        for i=1:LL-N*q
            if dmask[i] == false
                lb[ i ] = -1e-6
                ub[ i ] =  1e-6
                x0[ i ] =  0.0            
            else
                x0[ i ] += (1.0 - 2.0*rand())*1e-3
            end
        end        
    end    
    return lb,ub
end
