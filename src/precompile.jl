import Base:precompile

precompile(read_fasta, (String, Float64, Symbol, Bool))
precompile(compute_weights,(Matrix{Int8},Int,Symbol))
precompile(read_fasta_alignment,(AbstractString,Float64))
precompile(plmdca_asym,(Array{Float64,2},Vector{Float64}))
precompile(minimize_pl_asym,(PlmAlg,PlmVar))
precompile(pl_site_grad!,(Vector{Float64},Vector{Float64}))
precompile(fillvecene!,(Vector{Float64},Vector{Float64},Int, Vector{Int}, Int, Int))
precompile(logsumexp,(Vector{Float64},))
precompile(l2norm_asym,(Array{Float64,1}, PlmVar))
precompile(compute_score,(Array{Float64,2}, PlmVar, Int))
