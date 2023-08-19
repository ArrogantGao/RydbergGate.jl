module RydbergGate

# this package is to build a universal Hamiltonian via Rydberg atoms systems
# the package contains the following parts

using LinearAlgebra, GenericTensorNetworks, Graphs, Yao, LuxorGraphPlot, KrylovKit

# (a) an interface to accept the full Hamiltonian as input (using Yao.jl)

include("setup.jl")

# (b) a solver that can find the WMIS and theirs relation to all possible states

# (c) a solver used to calculate the off-diagional terms of the effective Hamiltonian (both analytical and numerical)

# (d) a comparation to direct results via Yao.jl


end
