# this script is used for setting up the system

using Yao.ConstGate: P1, P0
export RydbergArray

struct RydbergArray{TE, TΔ, Tδ, TΩ}
    N::Int # number of the vertices
    edges::Vector{Tuple{Int, Int}} # the edgess
    edges_binary::Vector{TE} # the edges in binary as Int64
    U::TΔ # blockade energy
    Δ::Vector{TΔ} # the excation energy
    δ::Vector{Tδ} # detuning on P1
    Ω::Vector{TΩ} # detuning on σₓ
end

function edges2binary(Type, edges::Vector{Tuple{Int, Int}})
    edges_binary = zeros(Type, length(edges))
    for n in 1:length(edges)
        i, j = edges[n]
        edges_binary[n] = edges_binary[n] | (Type(1) << i)
        edges_binary[n] = edges_binary[n] | (Type(1) << j)
    end
    return edges_binary
end

function type_choice(N::Int)
    if N < 8
        return UInt8
    elseif N < 16
        return UInt16
    elseif N < 32
        return UInt32
    elseif N < 64
        return UInt64
    elseif N < 128
        return UInt128
    else
        error("Too many atoms (more than 127)!")
    end
end

function RydbergArray(N::Int, edges::Vector{Tuple{Int, Int}}, U::TΔ, Δ::TΔ) where{TΔ <: Number}
    TE = type_choice(N)
    edges_binary = edges2binary(TE, edges)
    return RydbergArray{TE, TΔ, TΔ, TΔ}(N, edges, edges_binary, U, Δ .* ones(TΔ, N), zeros(TΔ, N), zeros(TΔ, N))
end

function RydbergArray(N::Int, edges::Vector{Tuple{Int, Int}}, U::TΔ, Δ::Vector{TΔ}) where{TΔ <: Number} 
    TE = type_choice(N)
    edges_binary = edges2binary(TE, edges)
    return RydbergArray{TE, TΔ, TΔ, TΔ}(N, edges, edges_binary, U, Δ, zeros(TΔ, N), zeros(TΔ, N))
end

function RydbergArray(N::Int, edges::Vector{Tuple{Int, Int}}, U::TΔ, Δ::Vector{TΔ}, δ::Vector{Tδ}) where{TΔ <: Number, Tδ <: Number} 
    TE = type_choice(N)
    edges_binary = edges2binary(TE, edges)
    return RydbergArray{TE, TΔ, Tδ, TΔ}(N, edges, edges_binary, U, Δ, δ, zeros(TΔ, N))
end

function RydbergArray(N::Int, edges::Vector{Tuple{Int, Int}}, U::TΔ, Δ::Vector{TΔ}, δ::Vector{Tδ}, Ω::Vector{TΩ}) where{TΔ <: Number, Tδ <: Number, TΩ <: Number} 
    TE = type_choice(N)
    edges_binary = edges2binary(TE, edges)
    return RydbergArray{TE, TΔ, Tδ, TΩ}(N, edges, edges_binary, U, Δ, δ, Ω)
end

# function to plot the built Rydberg Array

function plotRydbergArray(Rydberg_array::RydbergArray{TE, TΔ, Tδ, TΩ}) where{TE <: Unsigned, TΔ <: Number, Tδ <: Number, TΩ <: Number} 
    graph = SimpleGraph(Rydberg_array.N)
    for (i, j) in Rydberg_array.edges add_edge!(graph, i, j) end
    show_graph(graph)
    # nothing
end

# converter from RydbergArray to Hamiltonian in Yao.jl

function direct_map_Hamiltonian(Rydberg_array::RydbergArray{TE, TΔ, Tδ, TΩ}) where{TE <: Unsigned, TΔ <: Number, Tδ <: Number, TΩ <: Number} 
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