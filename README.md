PlmDCA
======
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pagnani.github.io/PlmDCA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pagnani.github.io/PlmDCA.jl/dev)
[![Build Status](https://github.com/pagnani/PlmDCA/workflows/CI/badge.svg)](https://github.com/pagnani/PlmDCA/actions) 
[![Coverage](https://codecov.io/gh/pagnani/PlmDCA/branch/master/graph/badge.svg)](https://codecov.io/gh/pagnani/PlmDCA)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Pseudo-likelihood maximization in [Julia](http://julialang.org). A
complete description of the algorithm can be found at
http://plmdca.csc.kth.se/. If you use this algorithm you should cite:

1. M. Ekeberg, C. Lovkvist, Y. Lan, M. Weigt, E. Aurell, Improved
   contact prediction in proteins: Using pseudolikelihood to infer Potts
   models, Phys. Rev. E 87, 012707 (2013)

2. M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
   maximization for direct-coupling analysis of protein structure from
   many homologous amino-acid sequences,
   [arXiv:1401.4832](http://arxiv.org/abs/1401.4832) (supplementary
   material)

The present software is a Julia implementation of above mentioned
papers, with no reference to the original MATLAB
[software](http://plmdca.csc.kth.se) implementation.

The code now requires at least Julia version 1.5 or later.

Install
-------

To install locally.
```
(v1.?) pkg> add https://github.com/pagnani/PlmDCA
```
The package is now registered and should soon appear in the General Registry. 

Overview
--------

The code internally uses [NLopt](https://github.com/JuliaOpt/NLopt.jl) which
provides a Julia interfaces to the free/open-source [NLopt
library](http://ab-initio.mit.edu/wiki/index.php/NLopt). 


Usage
-----
To load the code just type
```
julia> using PlmDCA
```

The functions in this package are written to maximize performance. Most
computationally-heavy functions can use multiple threads (start julia with 
the `-t` option or set the `JULIA_NUM_THREADS` environment variable). 
For more information on how set correctly the number of threads, please
refer to the online [Julia Documentation on Multi-Threading](https://docs.julialang.org/en/v1/manual/multi-threading/).

The software provides two main functions `plmdca(filename::String,
...)` and `plmdca_sym(filename::String,...)` (resp. the asymmetric and
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

After some experiments we found that the best compromise between
accuracy and speed is achieved by the Low Storage BFGS method
`:LD_LBFGS`, which is the default method in the code. The other
methods can be set changing the default optional argument
(e.g. `method=:LD_SLSQP`).

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
*  `X.score`:  a  vector of `Tuple{Int,Int,Float64}` containing the candidate contacts in descending score order (residue1, residue2 , score12).

Requirements
---

* The minimal julia version for using this code is 1.3 (package version <= v0.2.0)

* From package versions 0.3.0 on the minimal `julia` requirement is 1.5 (although
  the oldest version we test is v1.6)

