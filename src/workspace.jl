using  LinearAlgebra, SparseArrays, StaticArrays, BenchmarkTools

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

function Hadamard!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    H =  sparse([1/sqrt(2)  1/sqrt(2)
                 1/sqrt(2) -1/sqrt(2)]) 
    M = spzeros(2,2)
    for  i ∈ 0:a-1
        if (i == 0) && (i !== bitpos-1)
            copy!(M,sparse(I(2)))
        elseif (i == 0) && (i == bitpos-1)
            copy!(M,H)
        else  
            if (i == bitpos-1)
                M = kron(H,M)
            else
                M = kron(sparse(I(2)),M)
            end
        end
    end
   return #push!(C.gates, M)
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
   return #push!(C.gates, M)
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
    I2 = sparse(convert(Matrix{ComplexF64},I(2)))

    conket1 = sparse([0.0+0.0*im 0.0+0.0*im;
                      0.0+0.0*im 1.0+0.0*im])
    conket0 = sparse([1.0+0.0*im 0.0+0.0*im
                      0.0+0.0*im 0.0+0.0*im])

    for i ∈ 0:a-1

        if (i==0==control-1)
            M1 = conket1
        elseif (i==0==target-1)
            M1 = X
        else
            if (i==control-1)
                M1 = kron(conket1,M1)
            elseif (i==target-1)
                M1 = kron(X,M1)
            else
                M1 = kron(I2,M1)
            end
        end
    end

    for i ∈ 0:a-1

        if (i==0==control-1)
            M2 = conket0
        elseif (i==0!==control-1)
            M2 = I2
        else
            if (i==control-1)
                M2 = kron(conket0,M2)
            else
                M2 = kron(I2,M2)
            end
        end
    end
    return #push!(C.gates,(@view(M1)+@view(M2)))

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

function shots(C::Circuit, nshots::Int)
    a    
end