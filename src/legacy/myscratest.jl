function test(site)
    filename = "data/PF00014.fasta";
    max_gap_fraction = 0.9;
    theta = :auto;
    remove_dups = true;
    method = :LD_LBFGS;
    verbose = true;
    epsconv = 1e-5;
    maxit = 1000;
    lambdaJ = 0.01;
    lambdaH = 0.01;
    W, Z, N, M, q = PlmDCA.read_fasta(filename, max_gap_fraction, theta, remove_dups);
    plmvar = PlmDCA.PlmVar(N, M, q, lambdaJ, lambdaH, Z, W);
    plmalg = PlmDCA.PlmAlg(method, verbose, epsconv, maxit)
    LL = (plmvar.N - 1) * plmvar.q2 + plmvar.q
    x0 = zeros(Float64, LL);
    g  = similar(x0)
    Jmat = zeros(LL, plmvar.N)
    opt = Opt(plmalg.method, length(x0))
    ftol_abs!(opt, plmalg.epsconv)
    xtol_rel!(opt, plmalg.epsconv)
    xtol_abs!(opt, plmalg.epsconv)
    ftol_rel!(opt, plmalg.epsconv)
    maxeval!(opt, plmalg.maxit)
    #PlmDCA.min_objective!(opt, (x, g) -> PlmDCA.optimfunwrapper(x, g, site, plmvar))
    PlmDCA.min_objective!(opt, (x, g) -> PlmDCA.pl_site_grad!(x,g,site,plmvar))
    #(minf, minx, ret) = optimize(opt, x0)
    @code_warntype PlmDCA.pl_site_grad!(x0,g,site,plmvar)
end