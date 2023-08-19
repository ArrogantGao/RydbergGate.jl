# this script is used for setting up the system

using Yao.ConstGate: P1, P0
export RydbergArray

struct RydbergArray{TΔ, Tδ, TΩ}
    N::Int # number of the vertices
    edges::Vector{Tuple{Int, Int}} # the edgess
    U::TΔ # blockade energy
    Δ::Vector{TΔ} # the excation energy
    δ::Vector{Tδ} # detuning on P1
    Ω::Vector{TΩ} # detuning on σₓ
end

RydbergArray(N::Int, edges::Vector{Tuple{Int, Int}}, U::TΔ, Δ::TΔ) where{TΔ <: Number} = RydbergArray{TΔ, TΔ, TΔ}(N, edges, U, Δ .* ones(TΔ, N), zeros(TΔ, N), zeros(TΔ, N))
RydbergArray(N::Int, edges::Vector{Tuple{Int, Int}}, U::TΔ, Δ::Vector{TΔ}) where{TΔ <: Number} = RydbergArray{TΔ, TΔ, TΔ}(N, edges, U, Δ, zeros(TΔ, N), zeros(TΔ, N))
RydbergArray(N::Int, edges::Vector{Tuple{Int, Int}}, U::TΔ, Δ::Vector{TΔ}, δ::Vector{Tδ}) where{TΔ <: Number, Tδ <: Number} = RydbergArray{TΔ, Tδ, TΔ}(N, edges, U, Δ, δ, zeros(TΔ, N))
RydbergArray(N::Int, edges::Vector{Tuple{Int, Int}}, U::TΔ, Δ::Vector{TΔ}, δ::Vector{Tδ}, Ω::Vector{TΩ}) where{TΔ <: Number, Tδ <: Number, TΩ <: Number} = RydbergArray{TΔ, Tδ, TΩ}(N, edges, U, Δ, δ, Ω)

# function to plot the built Rydberg Array

function plotRydbergArray(Rydberg_array::RydbergArray{TΔ, Tδ, TΩ}) where{TΔ <: Number, Tδ <: Number, TΩ <: Number} 
    graph = SimpleGraph(Rydberg_array.N)
    for (i, j) in Rydberg_array.edges add_edge!(graph, i, j) end
    show_graph(graph)
    # nothing
end

# converter from RydbergArray to Hamiltonian in Yao.jl

function direct_map_Hamiltonian(Rydberg_array::RydbergArray{TΔ, Tδ, TΩ}) where{TΔ <: Number, Tδ <: Number, TΩ <: Number} 
    vertices = 1:Rydberg_array.N
    h = sum([Rydberg_array.U * kron(5, i=>P1, j=>P1) for (i, j) in Rydberg_array.edges]) + 
        sum([( - Rydberg_array.Δ[i] + Rydberg_array.δ[i]) * put(Rydberg_array.N, i=>P1) for i in vertices]) +
        sum([Rydberg_array.Ω[i] * put(Rydberg_array.N, i=>X) for i in vertices])
    return h
end

# this function is used to generate the eigen of the system directly
function direct_eigen_Hamiltonian(h::T) where{T <: Add}
    return eigen(Matrix(mat(h))).values
end

function direct_eigen_Hamiltonian(h::T, n::Int) where{T <: Add}
    return eigsolve(mat(h), randn(ComplexF64, 32), n, :SR)
end