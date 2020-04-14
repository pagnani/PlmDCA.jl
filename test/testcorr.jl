module TestCorr
using PlmDCA, Test

function corrtests_simple()
    W,Z,N,M,q = PlmDCA.ReadFasta("test.fasta",1.0, 0.0 , true)
    @test q == 2
    @test N == 3
    @test M == 2
    @test Z == [1 1; 2 2; 1 2]
    @test W == [0.5, 0.5]
end
corrtests_simple()
printstyled("All TestCorr passed!\n",color=:light_green,bold=true)
end # end module TestGauge
