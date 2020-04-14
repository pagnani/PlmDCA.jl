using PlmDCA,Test

tests = ["testdca","testcorr"]

printstyled("Running tests:\n", bold=true, color=:light_blue)

using Random
Random.seed!(345679)

res = map(tests) do t
    @eval module $(Symbol("Test_",t))
    using Random
    Random.seed!(345679)
    include($t * ".jl")
    end
    return
end
# print method ambiguities

ambiguities=Test.detect_ambiguities(PlmDCA)
printstyled("Potentially stale exports: ",bold=true, color=:light_blue)
isempty(ambiguities) ? (printstyled("\t none\n",bold=true, color=:light_green)) : (display(ambiguities)) 


