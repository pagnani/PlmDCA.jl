PlmDCA
======

Pseudo-likelihood maximization in [Julia](http://julialang.org). A
complete description of the algorithm can be found at
http://plmdca.csc.kth.se/. If you use this algorithm you should cite:

1. M. Ekeberg, C. Lovkvist, Y. Lan, M. Weigt, E. Aurell, Improved
   contact prediction in proteins: Using pseudolikelihoods to infer Potts
   models, Phys. Rev. E 87, 012707 (2013)

2. M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
   maximization for direct-coupling analysis of protein structure from
   many homologous amino-acid sequences,
   [arXiv:1401.4832](http://arxiv.org/abs/1401.4832) (supplementary
   material)

The present software is a Julia implementation of above mentioned
papers, with no reference to the original MATLAB
[software](http://plmdca.csc.kth.se) implementation.

The code now requires at least Julia version 0.6 or later.

Overview
--------

The code uses [NLopt](https://github.com/JuliaOpt/NLopt.jl) which
provides a Julia interfaces to the free/open-source [NLopt
library](http://ab-initio.mit.edu/wiki/index.php/NLopt). The program
(only in its asymmetric version so far, see below) can be run on
multiple cores previous ``addprocs(nprocs)`` where ``nprocs`` should
be some integer number `np` lower or equal to your (physical) number
of cores.

Install
-------

It requires the installation of:

* [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl). 
```
julia> Pkg.add("NLopt")
```
* [GaussDCA](https://github.com/carlobaldassi/GaussDCA.jl)
```
julia> Pkg.clone("https://github.com/carlobaldassi/GaussDCA.jl")
```
* [PlmDCA] The PlmDCA module itself can be added with
```
julia> Pkg.clone("https://github.com/pagnani/PlmDCA")
```
To install the package under version >= 0.7 please use the new package manager (it can be activated from the REPL using the key `]`)

```
(v0.7) pkg> add https://github.com/pagnani/PlmDCA #v07
```
We have not tested yet the software on Windows.

Usage
-----
To load the code just type
```
julia> using PlmDCA
```

The software provides two main functions `plmdca(filename::String,
...)` and `plmdcasym(filename::String,...)` (resp. the asymmetric and
symmetric coupling version of the algorithm). Empirically it turns out
that the asymmetric version is faster and more accurate. This function
take as input the name of a (possibly zipped) multiple sequence.

We also provide another function `mutualinfo(filename::String,...)` to
compute the mutual information score. 

There are a number of possible algorithmic strategies for the
optimization problem. As long as local gradient-based optimization is
concerned, this is a list of `:symbols` (associated to the different
methods): 
```
:LD_MMA, :LD_SLSQP, :LD_LBFGS, :LD_TNEWTON_PRECOND
:LD_TNEWTON_PRECOND_RESTART, :LD_TNEWTON, :LD_VAR2, :LD_VAR1
```

After some test we found that the best compromise between accuracy and
speed is achieved by the Low Storage BFGS method `:LD_LBFGS`, which is
the default method in the code. The other methods can be set changing
the default optional argument (e.g. `method=:LD_SLSQP`).

There are more optional arguments that can be set (to be documented...).

Output
======
The functions output a `type PlmOut` (say `X`) with 4 fields:


*  `X.Jtensor`: the coupling matrix `J[ri,rj,i,j]` a symmetrized
`q x q x N x N` array, where `N` is the number of residues in the
multiple sequence alignment, and `q` is the alphabet "size" (typically
21 for proteins).
*  `X.htensor`: the external field `h[r_i,i]` `q x N` array.
*  `X.pslike`: the pseudolikelihood
*  `X.score`: a `(Int, Int, Float64)` (in julia version 0.3) or a
`Tuple{Int,Int,Float64}` in (in julia version 0.4) vector of Tuples containing the candidate contact descending order (residue1, residue2 , score12).

Requirements
---

The minimal julia version for using this code is 0.6.

Todos 
----- 

A lot!
