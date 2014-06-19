function DecimateSym!(var::PlmVar, alg::PlmAlg, dec::DecVar{1})

    Jvec,pseudolike = MinimizePLSym(alg, var) 
    Nc2 = binomial(var.N,2)
    LL = Nc2 * var.q2
#    results = (Int,Float64,Float64,Float64)[]
    iteration = 1
    while true        
        SetAllCouplingToZero!(Jvec,  var, dec)
        nzeros = sum(!dec.dmask)
        fraczeros = nzeros / LL 
        println("# set to zero ", nzeros ," over ", LL ," equal to ", 100.0*fraczeros," % of the couplings. Converging ...")
        pseudolike, Jvec = MinimizePLSym(alg, var, dec, Jvec, x0=Jvec)
        println("# Converged!")
#        push!(results,(iteration,1.0-fraczeros, fmin,fminloc))
        fraczeros >= dec.fracmax && break
        iteration += 1
    end
    return Jvec, pseudolike
end

function SetAllCouplingToZero!(vecJ::Array{Float64,1}, var::PlmVar, dec::DecVar{1})

    N = var.N
    q  = var.q
    q2 = var.q2

    LL = length(vecJ)
    @assert LL == (binomial(N,2)*q2 + N*q)
    println("LL = $LL dmask = ", length(dec.dmask))
   
    vecval = (Int,Float64)[]
    sizehint(vecval,LL)
    for i=1:LL-N*q
        dec.dmask[i] && push!(vecval,(i, vecJ[i]*vecJ[i]))
    end
    
    
    lvec = length(vecval)
    idxperm = sortperm(Float64[ vecval[i][2] for i=1:lvec]) 
    
    ndec = int(dec.fracdec * (LL-N*q))
    for i=1:ndec
        site = vecval[idxperm[i]][1]
        dec.dmask[site] = false
    end
end




function SetCouplingToZero!(FN::Array{Float64,2},var::PlmVar, dec::DecVar{1})

    
    N = size(FN,1)
    
    q  = var.q
    q2 = var.q2

   
    vecval = (Int,Float64)[]
    counter = 0
    for i=1:(N-1), j=(i+1):N
        counter += 1
        dec.dmask[counter] && push!(vecval,(counter, FN[i,j])) 
    end
    
    
    lvec = length(vecval)
    idxperm = sortperm(Float64[ vecval[i][2] for i=1:lvec]) 
    
    ndec = int(dec.fracdec * counter)
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

    function f(x::Vector, g::Vector)
        eltime = @elapsed begin
        g === nothing && (g = zeros(Float64, length(x)))
        a = PLsiteAndGradSym!(x, g, var)
        end
        println("time spent in iteration = $eltime . PL = $a")
        return a
    end
        
    opt = Opt(alg.method, length(x0))

    ftol_abs!(opt, alg.epsconv)

    lb,ub = ComputeConstraintsSym!(dec, var.q, var.N, x0)

    lower_bounds!(opt, lb)

    upper_bounds!(opt, ub)

    maxeval!(opt, alg.maxit)

    min_objective!(opt, f)
    elapstime = @elapsed  (minf, minx, ret) = optimize(opt, x0)

    alg.verbose && @printf("pl = %.4f\t time = %.4f\n",  minf, elapstime)
    return minf,minx
end 

function ComputeConstraintsSym!(decvar::DecVar{1}, q::Int, N::Int, x0::Array{Float64,1} )
    
    q2 = q*q    
    dmask = decvar.dmask 
    Nc2 = binomial(N,2)
    lb = -Inf * ones(Float64, Nc2*q2 + N*q)
    ub =  Inf * ones(Float64, Nc2*q2 + N*q)
    LL = length(x0)

    for i=1:LL-N*q
        if dmask[i] == false
            lb[ i ] = -1.0e-6
            ub[ i ] =  1.0e-6
            x0[ i ] =  0.0            
        else
            x0[ i ] += (1.0 - 2.0*rand())*1e-3
        end
    end

    return lb,ub
end





