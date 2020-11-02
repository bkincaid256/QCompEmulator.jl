using  LinearAlgebra, SparseArrays, StaticArrays, BenchmarkTools

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

function Hadamard!(ψ::Qubit)
    #Hadamard gate for a single qubit    
    H = @SMatrix [1/sqrt(2)  1/sqrt(2)
                  1/sqrt(2) -1/sqrt(2)] 
    ψ.vec .= H*ψ.vec
    return ψ
end

function H!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    H =  sparse([1/sqrt(2)+0.0*im  1/sqrt(2)+0.0*im
                 1/sqrt(2)+0.0*im -1/sqrt(2)+0.0*im]) 
    M = spzeros(2,2)
    I2 = sparse(I*(1.0+0.0*im),2,2)
    for  i ∈ 1:a
        if (i == 1) && (i !== bitpos)
            copy!(M,I2)
        elseif (i == 1) && (i == bitpos)
            copy!(M,H)
        else  
            if (i == bitpos)
                M = kron(H,M)
            else
                M = kron(I2,M)
            end
        end
    end
   return push!(C.gates, M)
end

function H!(C::Circuit, set::Vector{Int})
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    for i in set
        H!(C,i)
    end
end

function H!(C::Circuit, set::UnitRange)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    for i in set
        H!(C,i)
    end
end

function Xgate!(ψ::Qubit)
    #X Pauli gate for a single qubit
    X =  @SMatrix [0.0 1.0
                   1.0 0.0]
    ψ.vec .= X*ψ.vec
    return ψ
end

function Xgate!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    X = sparse([0.0+0.0*im  1.0+0.0*im
                1.0+0.0*im  0.0+0.0*im]) 
    M = similar(X)
    I2 = sparse(I*(1.0+0.0*im),2,2)
    for  i ∈ 1:a
        if (i == 1) && (i !== bitpos)
            copy!(M,I2)
        elseif (i == 1) && (i == bitpos)
            copy!(M,X)
        else  
            if (i == bitpos)
                M = kron(X,M)
            else
                M = kron(I2,M)
            end
        end
    end
   return push!(C.gates, M)
end

function CX!(C::Circuit,control::Int,target::Int)

    X = sparse([0.0+0.0*im 1.0+0.0*im;
                1.0+0.0*im 0.0+0.0*im])

    a = length(C.qreg)

    if (target==control)
        error("Target qubit and control qubit must be Different!")
    elseif (target>a) || (control>a)
        error("The target and control qubit must be within the size of the register!")
    end

    M1 = spzeros(ComplexF64,2,2)
    M2 = copy(M1)
    I2 = sparse(I*(1.0+0.0*im),2,2)

    conket1 = sparse([0.0+0.0*im 0.0+0.0*im;
                      0.0+0.0*im 1.0+0.0*im])
    conket0 = sparse([1.0+0.0*im 0.0+0.0*im
                      0.0+0.0*im 0.0+0.0*im])

    for i ∈ 1:a

        if (i==1)
            if (i==control)
                M1 = conket1
            elseif (i==target)
                M1 = X
            else
                M1 = I2
            end
        else
            if (i==control)
                M1 = kron(conket1,M1)
            elseif (i==target)
                M1 = kron(X,M1)
            else
                M1 = kron(I2,M1)
            end
        end
    end

    for i ∈ 1:a

        if (i==1==control)
            M2 = conket0
        elseif (i==1!==control)
            M2 = I2
        else
            if (i==control)
                M2 = kron(conket0,M2)
            else
                M2 = kron(I2,M2)
            end
        end
    end
    return push!(C.gates,M1+M2)

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

function shots(C::Circuit, nshots::Int)
    a    
end