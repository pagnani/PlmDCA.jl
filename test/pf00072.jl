
using Pkg
Pkg.activate(".")
#using AnnealDCA

using GaussDCA
using PlmDCA

#using Test, Random, Statistics
#using StatsBase
#using PottsGauge

##########

filename="/home/guido/PROJECTS/pfam_pdb/data_pfam/PF00072/PF00072.fasta"

##########

Z = GaussDCA.read_fasta_alignment(filename, 0.8)
Z, _ = GaussDCA.remove_duplicate_seqs(Z)

N, M = size(Z)
q = round(Int,maximum(Z))
    
W , Meff = GaussDCA.compute_weights(Z,q,0.1)
W=W./sum(W)
#rmul!(W, 1.0/Meff)
Z=round.(Int,Z)



##############################
#Asym

#@time  resAplm = AnnealDCA.plmdca(Z,W);

time = @elapsed  resplm = PlmDCA.plmdca(Z,W);
#new -> 2088.275217671 
#old -> 2158.115226234



using JLD2

JLD2.@save "test/pf00072.jld2" resplm time

#= #############################
#Sym

@time  resAplm_s = AnnealDCA.plmdcasym(Z,W);

@time  resplm_s = PlmDCA.plmdca_sym(Z,W);

Js_old,hs_old=PottsGauge.gauge(resplm_s.Jtensor,resplm_s.htensor,PottsGauge.ZeroSumGauge())

println(sum(abs.(hs_old.-resAplm_s.htensor))/length(resAplm_s.htensor))

println(sum(abs.(Js_old.-resAplm_s.Jtensor))/length(resAplm_s.Jtensor))

#############################
 =#
#= 
using PyPlot

figure()

subplot(3,1,1)
plot(hs_old,resAplm_s.htensor,".")
subplot(3,1,2)
plot(hs_old,h_old,".")
subplot(3,1,3)
plot(h_old,resAplm_s.htensor,".")

gcf() =#