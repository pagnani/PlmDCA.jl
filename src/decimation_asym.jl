function DecimateAsym!(var::PlmVar, alg::PlmAlg, dec::DecVar{2})

    Jmat,pseudolike = MinimizePLAsym(alg,var) #initial minimization
    score , Jtensor = ComputeScore(Jmat, var)  #initial score

    LL,N = size(dec.dmask)
    allcoup = N*LL


    while true        
        SetCouplingToZeroAsym!(Jmat,  var, dec)
        nzeros = sum(!dec.dmask)
        fraczeros = nzeros/allcoup
        Jmat,pseudolike = MinimizePLAsym(alg, var, dec, Jmat, x0=Jmat)
        score, Jtensor = ComputeScore(Jmat, var)  #initial score
        fraczeros >= dec.fracmax && break
    end
    return Jmat,pseudolike
end

function SetCouplingToZeroAsym!(vecJ::Array{Float64,2}, var::PlmVar, dec::DecVar{2})

    LL,N = size(dec.dmask)
    q  = var.q
    q2 = var.q2
   
    vecval = (Int,Int,Float64)[]
    sizehint(vecval, N*LL)
    for j=1:N,i=1:LL
        dec.dmask[i,j] && push!(vecval,(i,j,vecJ[i,j]*vecJ[i,j])) # sorting absval
    end
    
    lvec = length(vecval)
    idxperm = sortperm(Float64[ vecval[i][3] for i=1:lvec]) 
    ndec = int(dec.fracdec * N*(N-1)*q2 )

    for i=1:ndec
        sitei = vecval[idxperm[i]][1]
        sitej = vecval[idxperm[i]][2]
        dec.dmask[sitei, sitej] = false
    end
end


function MinimizePLAsym(alg::PlmAlg, var::PlmVar, dec::DecVar{2}, Jmat::Array{Float64,2}; x0 = zeros(Float64,size(Jmat)))

    LL = size(x0,1)
    psvec = SharedArray(Float64, var.N)
    Jmatnew = @parallel hcat for site=1:var.N 
        function f(x::Vector, g::Vector)
            g === nothing  && (g = zeros(Float64, length(x)))
            return PLsiteAndGrad!(x, g, site,  var)            
        end
        vec0 = x0[:,site]
        lb,ub = ComputeConstraintsAsym!(dec, site, var.q, var.N, vec0)
        opt = Opt(alg.method, LL)
        ftol_abs!(opt, alg.epsconv)

        lower_bounds!(opt, lb)
        upper_bounds!(opt, ub)

        maxeval!(opt, alg.maxit)
        min_objective!(opt, f)
        elapstime = @elapsed  (minf, minx, ret) = optimize(opt, vec0)
        psvec[site] = minf
        alg.verbose && @printf("site = %d\t pl = %.4f\t time = %.4f\n", site, minf, elapstime)
        minx
    end 
    return Jmatnew, sum(psvec)
end 

function ComputeConstraintsAsym!(decvar::DecVar{2}, site::Int, q::Int, N::Int, x0::Array{Float64,1} )
    

    dmask = decvar.dmask 

    LL = length(x0)
    lb = -Inf * ones(Float64,LL)
    ub =  Inf * ones(Float64,LL)
    
    
    for i=1:size(decvar.dmask,1)
        if dmask[i,site] == false
            lb[ i ] = -1.0e-6
            ub[ i ] =  1.0e-6            
            x0[ i ] = 0.0 
        else
            x0[ i ] += (1.0 - 2.0*rand())*1e-3 # minimal pert to coupling 
        end
    end
    return lb,ub
end

