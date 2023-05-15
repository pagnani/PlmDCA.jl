using PrecompileTools

@setup_workload begin
    N = 10;
    M = 10
    Z = rand(1:21,N,M)
    W = rand(M)
    W ./= sum(W)
    @compile_workload begin
        redirect_stdout(devnull) do
            res1=plmdca_asym(Z,W)
            res2=plmdca_sym(Z,W)
            res3=mutualinfo(Int8.(Z))
        end
    end
end
