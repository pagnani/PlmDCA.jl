module PlmDCA
import Printf:@printf
import LinearAlgebra:norm, rmul!
import Statistics:mean
using NLopt
import DCAUtils: read_fasta_alignment, remove_duplicate_sequences, compute_weights, add_pseudocount, compute_weighted_frequencies
using LoopVectorization

export PlmOut, plmdca, plmdca_asym, plmdca_sym, plmdca_asym, mutualinfo

include("types.jl")
include("utils.jl")
include("plmdca_asym.jl")
include("plmdca_sym.jl")
include("mi.jl")
include("precompile.jl")
#include("legacy/decimation_sym.jl")
end
