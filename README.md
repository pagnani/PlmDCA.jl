plmdca.jl
========

Pseudo-likelihood maximization in [Julia](http://julialang.org). The
detail of the algorithm can be found at http://plmdca.csc.kth.se/

Overview
--------
The code uses [NLopt](https://github.com/JuliaOpt/NLopt.jl) which
provides a Julia interfaces to the free/open-source [NLopt
library](http://http://ab-initio.mit.edu/wiki/index.php/NLopt). The
program can be run on multiple cores previus ``addprocs(nprocs)``
where ``nprocs`` should be some integer number lower of equal your
(physical) number of cores.

Install
-------

It requires the installation of:
   
1. [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl). For Mac users with 10.9
   there are some known issues: try the following workaround:
```
julia> addprocs(nprocs)
julia> @everywhere dlopen("libstdc++",RTLD_GLOBAL)
julia> require("plmdca.jl")
```    
For Linuxes the installation of NLopt should be smooth.

2. [GaussDCA](https://github.com/carlobaldassi/GaussDCA.jl).

Not tested yet on Windows.

Usage
----
To load the code just type `require("plmdca.jl")`

The software provides one main function `plmdca(filename::String,
...)`. This function take as input the name of a (possibly zipped)
multiple sequence alignment in FASTA format, and returns 2 matrices:

* A matrix containing the coupling that minimize the Pseudolikelihood
function `J(col1, col2, pairij)` where `col1(2)` are the colors
(integers in the interval 1:21), and pairij is an integer in the
interval 1:N(N-1)/2, where N is the number of residues in the multiple
sequence alignment.

* A N x N matrix (N being the number of residues in the
  multiple-sequence alignment) with the scores. The higher the score
  the higher is the probability that residues are in contact.

Todos
-----
A lot. Make it a package.

