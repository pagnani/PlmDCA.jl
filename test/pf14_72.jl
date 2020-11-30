
using Pkg
Pkg.activate(".")
#using AnnealDCA

using GaussDCA
using PlmDCA
using BenchmarkTools
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

##############

@btime resplm = PlmDCA.plmdca(Z,W);

#############

filename="/home/guido/PROJECTS/pfam_pdb/data_pfam/PF00072/PF00072.fasta"

##########

Z = GaussDCA.read_fasta_alignment(filename, 0.8);
Z, _ = GaussDCA.remove_duplicate_seqs(Z);

N, M = size(Z)
q = round(Int,maximum(Z));
    
W , Meff = GaussDCA.compute_weights(Z,q,0.1);
W=W./sum(W);
#rmul!(W, 1.0/Meff)
Z=round.(Int,Z);

###########

@time resplm = PlmDCA.plmdca(Z,W);

