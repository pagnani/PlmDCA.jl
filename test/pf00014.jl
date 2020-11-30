
using Pkg
Pkg.activate(".")
#using AnnealDCA

using GaussDCA
using PlmDCA

#using Test, Random, Statistics
#using StatsBase
#using PottsGauge

##########

filename="/home/guido/PROJECTS/pfam_pdb/data_pfam/PF00014/PF00014.fasta"

##########

Z = GaussDCA.read_fasta_alignment(filename, 0.8)
Z, _ = GaussDCA.remove_duplicate_seqs(Z)

N, M = size(Z)
q = round(Int,maximum(Z))
    
W , Meff = GaussDCA.compute_weights(Z,q,0.1)

W=W./sum(W)
#rmul!(W, 1.0/Meff)
Z=round.(Int,Z)

@time resplm = PlmDCA.plmdca(Z,W);

#PF00014
#old plmDCA
#  23.376474 seconds (1.94 M allocations: 567.196 MiB, 0.14% gc time)
#new @avx
# 16.390855 seconds (10.66 M allocations: 1.945 GiB, 1.58% gc time)


#PF00072

#old plmDCA
#2158.115226234
#new @avx
#1212.079361911
