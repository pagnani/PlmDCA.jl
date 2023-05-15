```@meta
CurrentModule = PlmDCA
```

# PlmDCA
* Pseudo-likelihood maximization for protein.

## Package Features

- Learn model from multiple sequence alignment


See the [Index](@ref index) for the complete list of documented functions and types.

## Overview

Protein families are given in form of multiple sequence alignments (MSA) $D = (a^m_i |i = 1,\dots,L;\,m = 1,\dots,M)$ of $M$ proteins of aligned length $L$. The entries $a^m_i$ equal either one of the standard 20 amino acids, or the alignment gap $–$. In total, we have $q = 21$ possible different symbols in D. The aim of unsupervised generative modeling is to earn a statistical model $P(a_1,\dots,a_L)$ of (aligned) full-length sequences, which faithfully reflects the variability found in $D$: sequences belonging to the protein family of interest should have comparably high probabilities, unrelated sequences very small probabilities.
Here we propose a computationally efficient approach based on pseudo-likelihood maximization. 

We start from the exact decomposition $P_i(a_i| a_{\setminus i})$

Here, we use the following parametrization:

$P(a_i | a_1,\dots,a_{i-1}) = \frac{\exp \left\{ h_i(a_i) + \sum_{j=1}^{i-1} J_{i,j}(a_i,a_j)\right\} }{z_i(a_1,\dots,a_{i-1})}\,,$

where:

$z_i(a_1,\dots,a_{i-1})= \sum_{a=1}^{q} \exp \left\{ h_i(a) + \sum_{j=1}^{i-1} J_{i,j}(a,a_j)\right\} \,,$

is a the normalization factor. In machine learning, this
parametrization is known as soft-max regression, the generalization of logistic regression to multi-class labels.

# Usage

The typical pipeline to use the package is:

* Compute PlmDCA parameters from a multiple sequence alignment:

``` 
julia> res=plmdca(filefasta; kwds...)
```


# Multithreading

To fully exploit the the multicore parallel computation, julia should be invoked with

```
$ julia -t nthreads # put here nthreads equal to the number of cores you want to use
```

If you want to set permanently the number of threads to the desired value, you can either create a default environment variable `export JULIA_NUM_THREADS=24` in your `.bashrc`. More information [here](https://docs.julialang.org/en/v1.6/manual/multi-threading/)



## Load PlmDCA package 

The following cell loads the package `PlmDCA` (*Warning*: the first time it takes a while)

* The `mypkgdir` variable should be set to your `path/to/package` dir.

We will use the PF00014 protein family available in `data/` folder of the package/

```
mypkgdir = normpath(joinpath(pwd(),".."))
datadir=joinpath(mypkgdir,"data") # put here your path
using Pkg
Pkg.activate(mypkgdir)
using PlmDCA
```
## Learn the autoregressive parameters

As a preliminary step, we learn the field and the coupling parameters $h,J$ from the MSA. To do so we use the `plmdca` method that return the parameters stored in `res::PlmOut` in the cell below.

The keyword arguments for the `plmdca` method are (with their default value):

* `epsconv::Real=1.0e-5` (convergence parameter)

* `maxit::Int=1000` (maximum number of iteration - don't change)

* `verbose::Bool=true` (set to `false` to suppress printing on screen)

* `method::Symbol=:LD_LBFGS` (optimization method)


```
res=plmdca("data/PF14/PF00014_mgap6.fasta.gz", verbose=false, lambdaJ=0.02,lambdaH=0.001);
```

## 1. Contact Prediction

We can compute the epistatic score for residue-residue contact prediction. To do so, we can use the `epistatic_score` method. The epistatic score is computed on any target sequence of the MSA. Empirically, it turns out the the final score does not depend much on the choice of the target sequence. 

The output is contained in a `Vector` of `Tuple` ranked in descending score order. Each `Tuple` contains $i,j,s_{ij}$ where $s_{ij}$ is the DCA score of the residue pair $i,j$. The residue numbering is that of the MSA, and not of the unaligned full sequences.

```
target_sequence = 1
score=epistatic_score(arnet,arvar,target_sequence)
```


## [Methods Reference](@id index)
```@index
```

```@autodocs
Modules = [PlmDCA]
```
