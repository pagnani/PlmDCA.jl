plmdca.jl
========

Pseudo-likelihood maximization in [Julia](http://julialang.org). The
detail of the algorithm can be found at http://plmdca.csc.kth.se/. If
you use this algorithm you should cite:

1. M. Ekeberg, C. Lovkvist, Y. Lan, M. Weigt, E. Aurell, Improved
   contact prediction in proteins: Using pseudolikelihoods to infer Potts
   models, Phys. Rev. E 87, 012707 (2013)

2. M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
   maximization for direct-coupling analysis of protein structure from
   many homologous amino-acid sequences, arXiv:1401.4832 (supplementary
   material)

The present software is just a reimplementation of the original MATLAB
[software](http://plmdca.csc.kth.se) in Julia.

Overview
--------
The code uses [NLopt](https://github.com/JuliaOpt/NLopt.jl) which
provides a Julia interfaces to the free/open-source [NLopt
library](http://ab-initio.mit.edu/wiki/index.php/NLopt). The
program can be run on multiple cores previus ``addprocs(nprocs)``
where ``nprocs`` should be some integer number lower of equal your
(physical) number of cores.

Install
-------

It requires the installation of:
   
1. [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl). For Mac users
   with OS X 10.9 there are some known issues with NLopt. In case you
   hace trouble installing it read the discussion 
   [here](https://github.com/JuliaOpt/NLopt.jl/issues/6). For making
   all workers use NLopt, try the following workaround: 
```
julia> addprocs(nprocs)
julia> @everywhere dlopen("libstdc++",RTLD_GLOBAL) 
julia> require("plmdca.jl") 
``` 
   For Linuxes the installation of NLopt should be smooth.
   
2. [GaussDCA](https://github.com/carlobaldassi/GaussDCA.jl)/

We have not tested yet the software on Windows.

Usage
----
To load the code just type `require("plmdca.jl")`

The software provides one main function `plmdca(filename::String,
...)`. This function take as input the name of a (possibly zipped)
multiple sequence alignment in FASTA format, and returns two matrices:

* A `21 x 21 x N(N-1)/2` matrix containing the couplings that maximize
  the Pseudolikelihood function `J(ai, aj, pairij)` where `ai, aj` are
  the amino acid types (integers in the interval 1:21) of residue pair
  `i,j`, and `pairij` is an integer in the interval `1:N(N-1)/2`,
  where `N` is the number of residues in the multiple sequence
  alignment.

* A symmetric `N x N matrix` (N being the number of residues in the
  multiple-sequence alignment) with the scores. The higher the score
  the higher is the probability that residues are in contact.

There are a number of possible algorithmic strategies for the
optimization problem. As long as local gradient-based optimization is
concernedl, this is a list of `:symbols` (associated to the different
methods): 
```
:NLOPT_LD_MMA, :LD_SLSQP, :LD_LBFGS, :LD_TNEWTON_PRECOND
:LD_TNEWTON_PRECOND_RESTART, :LD_TNEWTON, :LD_VAR2, :NLOPT_LD_VAR1
```

After some test we found that the best compromise betwee accuracy and
speed is achieved by the Low Storage BFGS method `:LD_LBFGS`, which is
the default method in the code.

Todos 
----- 

A lot. Make it a package.

