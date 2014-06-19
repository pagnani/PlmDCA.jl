function plmdca(filename::String;
                decimation::Bool=false,
                fracmax::Float64 = 0.3,
                fracdec::Float64 = 0.1,
                remove_dups::Bool = true,
                max_gap_fraction::Real = 0.9, 
                theta = :auto, 
                lambdaJ::Real=0.005, 
                lambdaH::Real=0.01,
                gaugecol::Int=-1,
                epsconv::Real=1.0e-5,
                maxit::Int=1000,
                verbose::Bool=true,
                method::Symbol=:LD_LBFGS)

    W,Z,N,M,q = ReadFasta(filename,max_gap_fraction, theta, remove_dups)

    plmalg = PlmAlg(method,verbose, epsconv ,maxit)
    plmvar = PlmVar(N,M,q,q*q,gaugecol,lambdaJ,lambdaH,Z,W)


    if decimation  == false
        Jmat, pslike = MinimizePLAsym(plmalg,plmvar)
    else
        decvar = DecVar{2}(fracdec, fracmax, ones(Bool, (N-1)*q*q, N)) 
        Jmat, pslike = DecimateAsym!(plmvar, plmalg, decvar)        
    end
    score, FNAPC, Jtensor = ComputeScore(Jmat, plmvar)

    return output = PlmOut{Float64,4}(pslike, Jtensor, score)
end
    

function MinimizePLAsym(alg::PlmAlg, var::PlmVar)

    x0 = zeros(Float64,(var.N - 1) * var.q2 + var.q)
    vecps = SharedArray(Float64,var.N)

    Jmat = @parallel hcat for site=1:var.N #1:12
        function f(x::Vector, g::Vector)
            g === nothing && (g = zeros(Float64, length(x)))
            return PLsiteAndGrad!(x, g, site,  var)            
        end
        opt = Opt(alg.method, length(x0))
        ftol_abs!(opt, alg.epsconv)
        maxeval!(opt, alg.maxit)
        min_objective!(opt, f)
        elapstime = @elapsed  (minf, minx, ret) = optimize(opt, x0)
        alg.verbose && @printf("site = %d\t pl = %.4f\t time = %.4f\n", site, minf, elapstime)
        vecps[site] = minf
        minx
    end 
    return Jmat, sum(vecps)
end

function PLsiteAndGrad!(vecJ::Array{Float64,1},  grad::Array{Float64,1}, site::Int, plmvar::PlmVar)

    LL = length(vecJ)
    q2 = plmvar.q2
    q = plmvar.q
    gaugecol = plmvar.gaugecol
    N = plmvar.N
    M = plmvar.M
    Z = sdata(plmvar.Z)
    W = sdata(plmvar.W)


    for i=1:LL-q
        grad[i] = 2.0 * plmvar.lambdaJ * vecJ[i]
    end
    for i=(LL-q+1):LL
       grad[i] = 2.0 * plmvar.lambdaH * vecJ[i]
    end 

    vecene = zeros(Float64,q)
    expvecenesunorm = zeros(Float64,q)
    pseudolike = 0.0
 
    @inbounds begin 
        for a = 1:M       
            fillvecene!(vecene, vecJ,site,a, q, Z,N)        
            norm = sumexp(vecene)
            expvecenesunorm = exp(vecene .- log(norm))
            pseudolike -= W[a] * ( vecene[Z[site,a]] - log(norm) )
            offset = 0         
            for i = 1:site-1 
                @simd for s = 1:q
                    grad[ offset + s + q * ( Z[i,a] - 1 ) ] += W[a] *  expvecenesunorm[s]
                end
                grad[ offset + Z[site,a] + q * ( Z[i,a] - 1 ) ] -= W[a] 
                offset += q2 
            end
	    for i = site+1:N 
                @simd for s = 1:q
                    grad[ offset + s + q * ( Z[i,a] - 1 ) ] += W[a] *  expvecenesunorm[s]
                end
                grad[ offset + Z[site,a] + q * ( Z[i,a] - 1 ) ] -= W[a] 
                offset += q2 
            end

            @simd for s = 1:q 
                grad[ offset + s ] += W[a] *  expvecenesunorm[s] 
            end
	    grad[ offset + Z[site,a] ] -= W[a] 	
        end
    end

   

    if 1 <= gaugecol <= q         
        offset = 0;
        @inbounds begin 
            for i=1:N-1
                for s=1:q
                    grad[offset + gaugecol + q * (s - 1)  ] = 0.0; # Gauge!!! set gradJ[a,q] = 0
                    grad[offset + s + q * (gaugecol - 1) ] = 0.0; # Gauge!!! set gradJ[q,a] = 0
                end
                offset += q2
            end
            grad[offset + gaugecol] = 0.0 # Gauge!!! set gradH[q] = 0
        end
    end

    pseudolike += L2norm_asym(vecJ, plmvar)

    return pseudolike 
end

function fillvecene!(vecene::Array{Float64,1}, vecJ::Array{Float64,1},site::Int, a::Int, q::Int, sZ::DenseArray{Int,2},N::Int)

    q2 = q*q   
    Z = sdata(sZ)
    @inbounds begin
        for l = 1:q
            offset::Int = 0
            scra::Float64 = 0.0
            for i = 1:site-1 # Begin sum_i \neq site J
                scra += vecJ[offset + l + q * (Z[i,a]-1)] 
                offset += q2 
            end
            # skipping sum over residue site
    	    for i = site+1:N
                scra += vecJ[offset + l + q * (Z[i,a]-1)] 
                offset += q2 
            end # End sum_i \neq site J
            scra += vecJ[offset + l] # sum H 
            vecene[l] = scra
        end
    end
end


function L2norm_asym(vec::Array{Float64,1}, plmvar::PlmVar)
    q = plmvar.q    
    N = plmvar.N
    lambdaJ = plmvar.lambdaJ
    lambdaH = plmvar.lambdaH

    LL = length(vec)

    mysum1 = 0.0
    @inbounds @simd for i=1:(LL-q)
        mysum1 += vec[i] * vec[i]
    end
    mysum1 *= lambdaJ

    mysum2 = 0.0
    @inbounds @simd for i=(LL-q+1):LL
        mysum2 += vec[i] * vec[i]
    end
    mysum2 *= lambdaH

    return mysum1+mysum2
end
