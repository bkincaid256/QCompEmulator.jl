###########Importing Libraries and Defining Structs############
using  LinearAlgebra, SparseArrays, StaticArrays, BenchmarkTools, Base, Plots

mutable struct Qubit #type for holding qubit information
    bitpos::Int
    vec::Vector{ComplexF64}
end

mutable struct Circuit
    qreg::Array{Qubit,1}
    gates::Vector{SparseMatrixCSC{ComplexF64,Int64}}
end

function Circuit(nqubits::Int)
    G = Vector{SparseMatrixCSC{ComplexF64,Int64}}()
    return Circuit(makeregister(nqubits),G)
end

function makeregister(nqubits::Int)
    #function for generating a quantum register of qubits
    zerovec = [1.0+0.0*im, 0.0+0.0*im]
    return [Qubit(i, copy(zerovec)) for i ∈ 1:nqubits]
end

function Base.max(set::AbstractRange)
    if set[1]<set[end]
        return set[end]
    elseif set[1]>set[end]
        return set[1]
    else
        return set[1]
    end
end

function Base.min(set::AbstractRange)
    if set[1]<set[end]
        return set[1]
    elseif set[1]>set[end]
        return set[end]
    else
        return set[1]
    end
end

########## Running the Circuit ###########
function Run_circ(C::Circuit)

    maxiter = length(C.gates)
    maxbit = length(C.qreg)

    M = copy(C.gates[1])
    ket = copy(C.qreg[1].vec)
    for i ∈ 2:maxbit
    ket = kron(C.qreg[i].vec,ket)
    end

    for i ∈ 2:maxiter
        M .= C.gates[i]*M
    end

    return M*ket
end

function Run_circ(C::Circuit)

    maxiter = length(C.gates)
    maxbit = length(C.qreg)

    M = copy(C.gates[1])
    ket = copy(C.qreg[1].vec)
    for i ∈ 2:maxbit
    ket = kron(C.qreg[i].vec,ket)
    end

    for i ∈ 2:maxiter
        M .= C.gates[i]*M
    end

    return M*ket
end

function shots(C::Circuit,nshots::Int)
    
    finalket = Run_circ(C)
    ntrials = rand(nshots)
    ensemble = zeros(Float64,length(finalket))
    trials!(ensemble, finalket, ntrials)
end

function shots(C::Circuit, nqubit::Int,nshots::Int)
    
    finalket = Run_circ(C)
    ntrials = rand(nshots)
    ensemble = zeros(Int,length(finalket))
    trials!(ensemble, finalket, ntrials)
end

function trials!(ensemble::AbstractVector{Float64}, ket::AbstractVector{ComplexF64}, probs::AbstractVector{Float64})
    nshots = length(probs)

    @inbounds for j ∈ 1:nshots
        local p = probs[j]
        local sumold = 0.0
        local sumnew = 0.0 
        
        @inbounds for i ∈ 1:length(ket)
            sumnew = sumold + real(conj(ket[i])*ket[i])
            if (p<=sumnew)
                ensemble[i] += 1.0
                break
            end
            sumold = sumnew
        end
    end
    return ensemble
end

function trials!(ensemble::AbstractVector{Float64}, ket::AbstractVector{ComplexF64}, nqubit::Int, probs::AbstractVector{Float64})
    nshots = length(probs)

    @inbounds for j ∈ 1:nshots
        local p = probs[j]
        local sumold = 0.0
        local sumnew = 0.0 
        
        @inbounds for i ∈ 1:length(ket)
            sumnew = sumold + real(conj(ket[i])*ket[i])
            if (p<=sumnew)
                ensemble[i] += 1.0
                break
            end
            sumold = sumnew
        end
    end
    return ensemble
end

function graph(data::AbstractVector)
    pyplot()
    normalize!(data)
    max = length(data)-1
    nqubits = convert(Int,log(2,max+1))
    xs = [string(i,base=2,pad=nqubits) for i ∈ 0:max]

    bar(xs,data, leg=false, title="Simulation of Circuit", xlabel="qubit readout", ylabel="Probabilities", dpi = 300)
end