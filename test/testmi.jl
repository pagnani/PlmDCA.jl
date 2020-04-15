module TestMI
using Test, PlmDCA
MI = PlmDCA.mutualinfo("test.fasta",pseudocount=0)
@test MI == zero(MI)
printstyled("All TestMI passed!\n",color=:light_green,bold=true)
end
