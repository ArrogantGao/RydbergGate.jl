# in this script, the WMIS will be generated, with all possible states
# all possible states will be treated as vertex and connected by the fliping relation
# a DFS will be used to find all possible roution from one vertex to another

export RydbergArrayState, state_generator, AllStates

# the structure RydbergArrayState contains the information of a single state, and each state should be a independent set
struct RydbergArrayState{TE <: Unsigned, TΔ}
    state::TE
    energy::TΔ
end

function RydbergArrayState(Rydberg_array::RydbergArray{TE, TΔ, Tδ, TΩ}, state::TE) where{TΔ <: Number, Tδ <: Number, TΩ <: Number, TE <: Unsigned}

    for edge in Rydberg_array.edges_binary
        if count_ones(edge & state) == 2
            return nothing
        end
    end

    energy = - sum(((state >> i) & 1) * Rydberg_array.Δ[i] for i in 1:Rydberg_array.N)

    return RydbergArrayState{TE, TΔ}(state, energy)
end

# this function will give all possible combination states as a Vector
# N is the length of array and s is the expected size of the set
function state_generator(Rydberg_array::RydbergArray{TE, TΔ, Tδ, TΩ}, s::T) where{TΔ <: Number, Tδ <: Number, TΩ <: Number, TE <: Unsigned, T <: Integer}
    # max_size = binomial(N, s)
    states = Vector{TE}()
    state = zero(TE)
    N = Rydberg_array.N

    state_iterator!(Rydberg_array, N, 1, s, state, states)

    return states
end

function state_iterator!(Rydberg_array::RydbergArray{TE, TΔ, Tδ, TΩ}, N::T, n::T, s::T, state::TE, states::Vector{TE}) where{TΔ <: Number, Tδ <: Number, TΩ <: Number, T <: Integer, TE <: Unsigned}
    if s == 0
        blockade = false
        for edge in Rydberg_array.edges_binary
            if count_ones(edge & state) == 2
                blockade = true
                break
            end
        end
        if blockade == false
            push!(states, state)
        end
    else
        for i in n : N - s + 1
            state_i = copy(state)
            state_i = state_i | (one(TE) << i)
            state_iterator!(Rydberg_array, N, i + 1, s - 1, state_i, states)
        end
    end
    return nothing
end

# this structure is used to store all the possible states.
# states will the same set size will be stored in the same Vector
# states_size will store the number of possible states
struct AllStates{T}
    states::Vector{Vector{RydbergArrayState{T}}}
    sizes::Vector{Int}
end

function AllStates(Rydberg_array::RydbergArray{TΔ, Tδ, TΩ}) where{TΔ <: Number, Tδ <: Number, TΩ <: Number}
    all_Rydberg_states = Vector{Vector{RydbergArrayState{TΔ}}}()
    states_size = Vector{Int}()
    for n in 0:Rydberg_array.N
        all_Rydberg_states_n = Vector{RydbergArrayState{TΔ}}()
        states_n = state_generator(Rydberg_array.N, n)
        for state_n in states_n
            rydberg_state = RydbergArrayState(Rydberg_array, state_n)
            if rydberg_state !== nothing
                push!(all_Rydberg_states_n, rydberg_state)
            end
        end

        size_n = length(all_Rydberg_states_n)
        if size_n != 0
            push!(all_Rydberg_states, all_Rydberg_states_n)
            push!(states_size, size_n)
        else
            break
        end
    end

    return AllStates{TΔ}(all_Rydberg_states, states_size)
end