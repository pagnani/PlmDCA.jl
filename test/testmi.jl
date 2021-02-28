module TestMI
using Test, PlmDCA
# MI = PlmDCA.mutualinfo("test.fasta",pseudocount=0)
#@test MI == zero(MI)
@test true
#printstyled("All TestMI passed!\n",color=:light_green,bold=true)
printstyled("TestMI waiting that https://github.com/carlobaldassi/DCAUtils.jl/issues/2 is fixed\n",color=:light_green,bold=true)
end
