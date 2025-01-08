var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PlmDCA","category":"page"},{"location":"#PlmDCA","page":"Home","title":"PlmDCA","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pseudo-likelihood maximization for protein.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Learn model from multiple sequence alignment","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Index for the complete list of documented functions and types.","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Protein families are given in form of multiple sequence alignments (MSA) D = (a^m_i i = 1dotsLm = 1dotsM) of M proteins of aligned length L. The entries a^m_i equal either one of the standard 20 amino acids, or the alignment gap . In total, we have q = 21 possible different symbols in D. The aim of unsupervised generative modeling is to earn a statistical model P(a_1dotsa_L) of (aligned) full-length sequences, which faithfully reflects the variability found in D: sequences belonging to the protein family of interest should have comparably high probabilities, unrelated sequences very small probabilities. Here we propose a computationally efficient approach based on pseudo-likelihood maximization. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"We start from the exact decomposition P_i(a_i a_setminus i)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here, we use the following parametrization:","category":"page"},{"location":"","page":"Home","title":"Home","text":"P(a_i  a_1dotsa_i-1) = fracexp left h_i(a_i) + sum_j=1^i-1 J_ij(a_ia_j)right z_i(a_1dotsa_i-1)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where:","category":"page"},{"location":"","page":"Home","title":"Home","text":"z_i(a_1dotsa_i-1)= sum_a=1^q exp left h_i(a) + sum_j=1^i-1 J_ij(aa_j)right ","category":"page"},{"location":"","page":"Home","title":"Home","text":"is a the normalization factor. In machine learning, this parametrization is known as soft-max regression, the generalization of logistic regression to multi-class labels.","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The typical pipeline to use the package is:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Compute PlmDCA parameters from a multiple sequence alignment:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> res=plmdca(filefasta; kwds...)","category":"page"},{"location":"#Multithreading","page":"Home","title":"Multithreading","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To fully exploit the the multicore parallel computation, julia should be invoked with","category":"page"},{"location":"","page":"Home","title":"Home","text":"$ julia -t nthreads # put here nthreads equal to the number of cores you want to use","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you want to set permanently the number of threads to the desired value, you can either create a default environment variable export JULIA_NUM_THREADS=24 in your .bashrc. More information here","category":"page"},{"location":"#Load-PlmDCA-package","page":"Home","title":"Load PlmDCA package","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The following cell loads the package PlmDCA (Warning: the first time it takes a while)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The mypkgdir variable should be set to your path/to/package dir.","category":"page"},{"location":"","page":"Home","title":"Home","text":"We will use the PF00014 protein family available in data/ folder of the package/","category":"page"},{"location":"","page":"Home","title":"Home","text":"mypkgdir = normpath(joinpath(pwd(),\"..\"))\ndatadir=joinpath(mypkgdir,\"data\") # put here your path\nusing Pkg\nPkg.activate(mypkgdir)\nusing PlmDCA","category":"page"},{"location":"#Learn-the-autoregressive-parameters","page":"Home","title":"Learn the autoregressive parameters","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"As a preliminary step, we learn the field and the coupling parameters hJ from the MSA. To do so we use the plmdca method that return the parameters stored in res::PlmOut in the cell below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The keyword arguments for the plmdca method are (with their default value):","category":"page"},{"location":"","page":"Home","title":"Home","text":"epsconv::Real=1.0e-5 (convergence parameter)\nmaxit::Int=1000 (maximum number of iteration - don't change)\nverbose::Bool=true (set to false to suppress printing on screen)\nmethod::Symbol=:LD_LBFGS (optimization method)","category":"page"},{"location":"","page":"Home","title":"Home","text":"res=plmdca(\"data/PF14/PF00014_mgap6.fasta.gz\", verbose=false, lambdaJ=0.02,lambdaH=0.001);","category":"page"},{"location":"#.-Contact-Prediction","page":"Home","title":"1. Contact Prediction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We can compute the epistatic score for residue-residue contact prediction. To do so, we can use the epistatic_score method. The epistatic score is computed on any target sequence of the MSA. Empirically, it turns out the the final score does not depend much on the choice of the target sequence. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The output is contained in a Vector of Tuple ranked in descending score order. Each Tuple contains ijs_ij where s_ij is the DCA score of the residue pair ij. The residue numbering is that of the MSA, and not of the unaligned full sequences.","category":"page"},{"location":"","page":"Home","title":"Home","text":"target_sequence = 1\nscore=epistatic_score(arnet,arvar,target_sequence)","category":"page"},{"location":"#index","page":"Home","title":"Methods Reference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PlmDCA]","category":"page"},{"location":"#PlmDCA.plmdca-Tuple{String}","page":"Home","title":"PlmDCA.plmdca","text":"plmdca(filename::String; kwds...) -> PlmOut\n\nRun plmdca_asym on the fasta alignment in filename It returns a struct ::PlmOut containing four fields:\n\nJtensor: the 4 dimensional Array (q×q×L×L)  of the couplings in zero-sum gauge J[a,b,i,j]\nHtensor: the 2 dimensional Array (q×L) of the fields in zero-sum gauge h[a,i]\npslike: a vector of length L containining the log-pseudolikelihoods\nscore: a vector of `(i,j,val)::Tuple{Int,Int,Float64} containing the DCA score val relative to the i,j pair in descending order.\nOptional arguments:\n\nlambdaJ::Real=0.01 coupling L₂ regularization parameter (lagrange multiplier)\nlambdaH::Real=0.01 field L₂ regularization parameter (lagrange multiplier)\nepsconv::Real=1.0e-5 convergence value in minimzation\nmaxit::Int=1000 maximum number of iteration in minimization\nverbose::Bool=true set to false to stop printing convergence info on stdout\nmethod::Symbol=:LD_LBFGS optimization strategy see NLopt.jl for other options\n\nExample\n\njulia> res = plmdca_asym(\"data/pf00014short.fasta\",lambdaJ=0.01,lambdaH=0.01,epsconv=1e-12);\n\n\n\n\n\n","category":"method"},{"location":"#PlmDCA.plmdca_asym-Tuple{String}","page":"Home","title":"PlmDCA.plmdca_asym","text":"plmdca_asym(filename::String; kwds...)-> PlmOut\n\nRun plmdca_asym on the fasta alignment in filename It returns a struct ::PlmOut containing four fields:\n\nJtensor: the 4 dimensional Array (q×q×L×L)  of the couplings in zero-sum gauge J[a,b,i,j]\nHtensor: the 2 dimensional Array (q×L) of the fields in zero-sum gauge h[a,i]\npslike: a vector of length L containining the log-pseudolikelihoods\nscore: a vector of `(i,j,val)::Tuple{Int,Int,Float64} containing the DCA score val relative to the i,j pair in descending order.\nOptional arguments:\n\nlambdaJ::Real=0.01 coupling L₂ regularization parameter (lagrange multiplier)\nlambdaH::Real=0.01 field L₂ regularization parameter (lagrange multiplier)\nepsconv::Real=1.0e-5 convergence value in minimzation\nmaxit::Int=1000 maximum number of iteration in minimization\nverbose::Bool=true set to false to stop printing convergence info on stdout\nmethod::Symbol=:LD_LBFGS optimization strategy see NLopt.jl for other options\n\nExample\n\njulia> res = plmdca_asym(\"data/pf00014short.fasta\",lambdaJ=0.01,lambdaH=0.01,epsconv=1e-12);\n\n\n\n\n\n","category":"method"},{"location":"#PlmDCA.plmdca_asym-Union{Tuple{T}, Tuple{Matrix{T}, Vector{Float64}}} where T<:Integer","page":"Home","title":"PlmDCA.plmdca_asym","text":"plmdca_asym(Z::Array{T,2},W::Vector{Float64}; kwds...) -> ::PlmOut\n\nPseudolikelihood maximization on the M×N alignment Z (numerically encoded in  1,…,21), and the M-dimensional normalized weight vector W.\n\nIt returns a struct ::PlmOut containing four fields:\n\nJtensor: the 4 dimensional Array (q×q×L×L)  of the couplings in zero-sum gauge J[a,b,i,j]\nHtensor: the 2 dimensional Array (q×L) of the fields in zero-sum gauge h[a,i]\npslike: a vector of length L containining the log-pseudolikelihoods\nscore: a vector of `(i,j,val)::Tuple{Int,Int,Float64} containing the DCA score val relative to the i,j pair in descending order.\nOptional arguments:\n\nlambdaJ::Real=0.01 coupling L₂ regularization parameter (lagrange multiplier)\nlambdaH::Real=0.01 field L₂ regularization parameter (lagrange multiplier)\nepsconv::Real=1.0e-5 convergence value in minimzation\nmaxit::Int=1000 maximum number of iteration in minimization\nverbose::Bool=true set to false to stop printing convergence info on stdout\nmethod::Symbol=:LD_LBFGS optimization strategy see NLopt.jl for other options\n\nExample\n\njulia> res = plmDCA(Z,W,lambdaJ=0.01,lambdaH=0.01,epsconv=1e-12);\n\n\n\n\n\n","category":"method"}]
}