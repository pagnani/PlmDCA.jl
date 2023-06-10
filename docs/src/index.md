```@meta
CurrentModule = PlmDCA
```

# PlmDCA
* Pseudo-likelihood maximization for protein.

## Package Features

- Learn pseudo-likelihood model from multiple sequence alignment


See the [Index](@ref index) for the complete list of documented functions and types.

## Overview

Protein families are given in form of multiple sequence alignments (MSA) $D =
(a^m_i |i = 1,\dots,L;\,m = 1,\dots,M)$ of $M$ proteins of aligned length $L$.
The entries $a^m_i$ equal either one of the standard 20 amino acids, or the
alignment gap $–$. In total, we have $q = 21$ possible different symbols in D.
The aim of unsupervised generative modeling is to earn a statistical model
$P(a_1,\dots,a_L)$ of (aligned) full-length sequences, which faithfully reflects
the variability found in $D$: sequences belonging to the protein family of
interest should have comparably high probabilities, unrelated sequences very
small probabilities. Here we propose a computationally efficient approach based
on pseudo-likelihood maximization. 

We start from the exact decomposition $P_i(a_i| a_{\setminus i})$ where
$a_{\setminus i} := \{a_1,\dots,a_{i-1},a_{i+1},\dots,a_i\}$, i.e. all residues of the sequence
of amino acids but the one relative to residue $i$. 

Here, we use the following parametrization:

$P(a_i |a_{\setminus i}) = \frac{\exp \left\{ h_i(a_i) + \sum_{j\neq i}
J_{i,j}(a_i,a_j)\right\} }{z_i(a_{\setminus i})}\,,$

where:

$z_i(a_{\setminus i})= \sum_{a=1}^{q} \exp \left\{ h_i(a) + \sum_{j=1\neq i} J_{i,j}(a,a_j)\right\} \,,$

is the normalization factor. The pseudo-likelihood maximization strategy, aims
at finding the value of the $J, h$ parameters, that maximize the log-likelihood:

${\cal L} = \frac{1}{M} \sum_{m=1}^M  \left( h_i(a^m_i) + \sum_{j\neq i}
J_{ij}(a_i^m,a_j^m) - \log z_i(a^m_{\setminus i})\right)\,.$

In machine learning, this parametrization is known as the {\em soft-max
regression}, a generalization of logistic regression to multi-class labels.

## Usage

The typical pipeline to use the package is:

* Compute PlmDCA parameters from a multiple sequence alignment:

``` 
julia> res=plmdca(filefasta; kwds...)
```

## Multithreading

To fully exploit the the multicore parallel computation, julia should be invoked with

```
$ julia -t nthreads # put here nthreads equal to the number of cores you want to use
```

If you want to set permanently the number of threads to the desired value, you can either create a default environment variable `export JULIA_NUM_THREADS=24` in your `.bashrc`. More information [here](https://docs.julialang.org/en/v1.6/manual/multi-threading/)

## Load PlmDCA package 

The package is on the General Registry. It can be installed from the package
manager by
```
pkg> add PlmDCA
```
and 
```
julia> using PlmDCA
```

## Learn the parameters

There are two different learning strategies: 

* The asymmetric one invoked by the `plmdca_asym` method (the `plmdca` method
  points to the asymmetric strategy)
* The symmetric strategy, invoked by the `plmdca_sym`. This method is slower and
  typically less accurate.

Both methods return the parameters  `res::PlmOut`, a `struct` containing:

1. `Jtensor`: the 4 dimensional Array (q×q×L×L)  of the couplings in zero-sum gauge J[a,b,i,j]
2. `Htensor`: the 2 dimensional Array (q×L) of the fields in zero-sum gauge h[a,i]
3. `pslike`: a vector of length `L` containing the log-pseudolikelihoods
4. `score`: a vector of ``(i,j,val)::Tuple{Int,Int,Float64}` containing the DCA
   score `val` relative to the `i,j` pair in descending order.

For both methods, the keyword arguments (with their default value) are:

* `epsconv::Real=1.0e-5` (convergence parameter)
* `maxit::Int=1000` (maximum number of iteration - don't change)
* `verbose::Bool=true` (set to `false` to suppress printing on screen)
* `method::Symbol=:LD_LBFGS` (optimization method)

```
res=plmdca("data/PF14/PF00014_mgap6.fasta.gz", verbose=false, lambdaJ=0.02,lambdaH=0.001);
```

## Contact Prediction

Contact prediction is contained in the the output of `plmdca`. `::Output`
contains a `score` field which is a  `Vector` of `Tuple` ranked in descending score order. Each `Tuple`
contains $i,j,s_{ij}$ where $s_{ij}$ is the DCA score of the residue pair $i,j$.
The residue numbering is that of the MSA, and not of the unaligned full sequences.

## [Methods Reference](@id index)
```@index
```

```@autodocs
Modules = [PlmDCA]
```
