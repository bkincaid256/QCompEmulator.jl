using  LinearAlgebra, StaticArrays, BenchmarkTools

mutable struct Qubit #type for holding qubit information
    bitpos::Int
    vec::Vector{ComplexF64}
end

mutable struct Circuit
    qreg::Array{Qubit,1}
    gates::Array{Array{Complex{Float64},2},1}
end

function makeregister(nqubits::Int)
    #function for generating a quantum register of qubits
    zerovec = [1.0+0.0*im, 0.0+0.0*im]
    return [Qubit(i, copy(zerovec)) for i ∈ 1:nqubits]
end

function Hadamard!(ψ::Qubit)
    #Hadamard gate for a single qubit    
    H = @SMatrix [1/sqrt(2)  1/sqrt(2)
                  1/sqrt(2) -1/sqrt(2)] 
    ψ.vec .= H*ψ.vec
    return ψ
end

function Hadamard!(C::Circuit, bitpos::Integer)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    H = @SMatrix [1/sqrt(2)  1/sqrt(2)
                  1/sqrt(2) -1/sqrt(2)] 
    M = Array{Float64}(undef,2,2)
    for  i ∈ 0:a-1
        if (i == 0) && (i !== bitpos-1)
            copy!(M,I(2))
        elseif (i == 0) && (i == bitpos-1)
            copy!(M,H)
        else  
            if (i == bitpos-1)
                M = kron(H,M)
            else
                M = kron(I(2),M)
            end
        end
    end
   return push!(C.gates, M)
end


function Xgate!(ψ::Qubit)
    #X Pauli gate for a single qubit
    X =  @SMatrix [0.0 1.0
                   1.0 0.0]
    ψ.vec .= X*ψ.vec
    return ψ
end

function Xgate!(C::Circuit, bitpos::Integer)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    X = @SMatrix [0.0  1.0
                  1.0  0.0] 
    M = Array{Float64}(undef,2,2)
    for  i ∈ 0:a-1
        if (i == 0) && (i !== bitpos-1)
            copy!(M,I(2))
        elseif (i == 0) && (i == bitpos-1)
            copy!(M,X)
        else  
            if (i == bitpos-1)
                M = kron(X,M)
            else
                M = kron(I(2),M)
            end
        end
    end
   return push!(C.gates, M)
end

function ControlX!(ψₜ::Qubit,ψc::Qubit)
    
    cX = @SMatrix [1.0 0.0 0.0 0.0
                   0.0 0.0 0.0 1.0
                   0.0 0.0 1.0 0.0
                   0.0 1.0 0.0 0.0]
    ket = kron(ψₜ.vec, ψc.vec)
    ket .= cX*ket     
    return ket
end

function ControlX!(C::Circuit,control::Integer,target::Integer)

    if (target==control)
        error("Target qubit and control qubit must be Different!")
    end

    if target > control
        cX = @SMatrix [1.0 0.0 0.0 0.0
                       0.0 0.0 0.0 1.0
                       0.0 0.0 1.0 0.0
                       0.0 1.0 0.0 0.0]

    else
        cX = @SMatrix [1.0 0.0 0.0 0.0
                       0.0 1.0 0.0 0.0
                       0.0 0.0 0.0 1.0
                       0.0 0.0 1.0 0.0]

    end
        a = length(C.qreg)
        M = Array{Float64}(undef,2,2)
        G = Matrix{Float64}[]
        if abs(target-control) == 1    

            bit = min(target,control)
            for i ∈ 0:a-2
                if (i == 0) && (bit-1 == 0)
                    M = kron(I(1),cX)
                elseif (i == 0) && (bit-1 !== 0)
                    M = kron(I(1),I(2))
                else
                    if (i == bit-1)
                        M = kron(cX,M)
                    else
                        M = kron(I(2),M)
                    end
                end
            end
        else
            ControlX!(C,target-1,control)
            ControlX!(C,target,control+1)
        end
    return push!(C.gates, M)
end

function Run_circ(C::Circuit)

    maxiter = length(C.gates::Array{Array{Complex{Float64},2},1})
    maxbit = length(C.qreg)

    M = convert(Matrix{ComplexF64},I(size(C.gates[begin])[1]))
    ket = copy(C.qreg[1].vec)
    for i ∈ 2:maxbit
    ket = kron(C.qreg[i].vec,ket)
    end

    for i ∈ 1:maxiter
        M .= C.gates[i]*M
    end

    return M*ket
end