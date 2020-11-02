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
######Single qubit gates ###########

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

function H!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        H!(C,i)
    end
end

function Ugate!(C::Circuit, U::AbstractArray{ComplexF64,2}, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    U = sparse(U)
    M = similar(U)
    I2 = sparse(I*(1.0+0.0*im),2,2)
    for  i ∈ 1:a
        if (i == 1) && (i !== bitpos)
            copy!(M,I2)
        elseif (i == 1) && (i == bitpos)
            copy!(M,U)
        else  
            if (i == bitpos)
                M = kron(U,M)
            else
                M = kron(I2,M)
            end
        end
    end
   return push!(C.gates, M)
end

function Ugate!(C::Circuit, U::AbstractArray{ComplexF64,2}, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        Ugate!(C,U,i)
    end
end

function T!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    T = sparse([1.0+0.0*im  0.0+0.0*im
                0.0+0.0*im  exp(im*π/4)]) 
    M = similar(T)
    I2 = sparse(I*(1.0+0.0*im),2,2)
    for  i ∈ 1:a
        if (i == 1) && (i !== bitpos)
            copy!(M,I2)
        elseif (i == 1) && (i == bitpos)
            copy!(M,T)
        else  
            if (i == bitpos)
                M = kron(T,M)
            else
                M = kron(I2,M)
            end
        end
    end
   return push!(C.gates, M)
end

function T!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        T!(C,i)
    end
end

######## Pauli Gates ##########
function X!(C::Circuit, bitpos::Int)
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

function X!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        X!(C,i)
    end
end

function Y!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    Y = sparse([0.0+0.0*im  0.0-1.0*im
                0.0+1.0*im  0.0+0.0*im]) 
    M = similar(Y)
    I2 = sparse(I*(1.0+0.0*im),2,2)
    for  i ∈ 1:a
        if (i == 1) && (i !== bitpos)
            copy!(M,I2)
        elseif (i == 1) && (i == bitpos)
            copy!(M,Y)
        else  
            if (i == bitpos)
                M = kron(Y,M)
            else
                M = kron(I2,M)
            end
        end
    end
   return push!(C.gates, M)
end

function Y!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        Y!(C,i)
    end
end

function Z!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    Z = sparse([1.0+0.0*im  0.0+0.0*im
                0.0+0.0*im -1.0+0.0*im]) 
    M = similar(Z)
    I2 = sparse(I*(1.0+0.0*im),2,2)
    for  i ∈ 1:a
        if (i == 1) && (i !== bitpos)
            copy!(M,I2)
        elseif (i == 1) && (i == bitpos)
            copy!(M,Z)
        else  
            if (i == bitpos)
                M = kron(Z,M)
            else
                M = kron(I2,M)
            end
        end
    end
   return push!(C.gates, M)
end

function Z!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        Z!(C,i)
    end
end

############# Controlled Gates and 2 qubit Gates#############
function CX!(C::Circuit,control::Int,target::Int)
    #Implentation of CNOT gate using Kronnecker products.
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

function Swap!(C::Circuit,b1::Int,b2::Int)
    #quick and dirty implmentation of Swap Gate
    CX!(C,b2,b1)
    CX!(C,b1,b2)
    CX!(C,b2,b1)
end

function CZ!(C::Circuit,control::Int,target::Int)
    #Implentation of CZ gate using Kronnecker products.
    Z = sparse([1.0+0.0*im  0.0+0.0*im;
                0.0+0.0*im -1.0+0.0*im])

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
                M1 = Z
            else
                M1 = I2
            end
        else
            if (i==control)
                M1 = kron(conket1,M1)
            elseif (i==target)
                M1 = kron(Z,M1)
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

function CCX!(C::Circuit,c1::Int,c2::Int,target::Int)
    #Implementation of the Toffoli gate using Kronnecker products. 
    X = sparse([0.0+0.0*im 1.0+0.0*im;
                1.0+0.0*im 0.0+0.0*im])

    a = length(C.qreg)

    if (target==c1) || (target==c2)
        error("Target qubit and control qubits must be Different!")
    elseif (c1==c2)
        error("The control qubits must be different from one another!")
    elseif ((target>a) || (c1>a) || (c2>a)) || ((target<=0) || (c1<=0) ||(c2<=0) )
        error("The target and control qubits must be within the size of the register!")
    end

    M1 = spzeros(ComplexF64,2,2)
    M2 = copy(M1)
    M3 = copy(M1)
    M4 = copy(M1)
    I2 = sparse(I*(1.0+0.0*im),2,2)

    conket1 = sparse([0.0+0.0*im 0.0+0.0*im;
                      0.0+0.0*im 1.0+0.0*im])
    conket0 = sparse([1.0+0.0*im 0.0+0.0*im
                      0.0+0.0*im 0.0+0.0*im])

    for i ∈ 1:a

        if (i==1)
            if (i==c1)||(i==c2)
                M1 = conket1
            elseif (i==target)
                M1 = X
            else
                M1 = I2
            end
        else
            if (i==c1)||(i==c2)
                M1 = kron(conket1,M1)
            elseif (i==target)
                M1 = kron(X,M1)
            else
                M1 = kron(I2,M1)
            end
        end
    end

    for i ∈ 1:a

        if (i==1)
            if (i==c1)||(i==c2)
                M2 = conket0
            else
                M2 = I2
            end
        else
            if (i==c1)||(i==c2)
                M2 = kron(conket0,M2)
            else
                M2 = kron(I2,M2)
            end
        end
    end

    for i ∈ 1:a

        if (i==1)
            if (i==c1)
                M3 = conket0
            elseif (i==c2)
                M3 = conket1
            else
                M3 = I2
            end
        else
            if (i==c1)
                M3 = kron(conket0,M3)
            elseif (i==c2)
                M3 = kron(conket1,M3)
            else
                M3 = kron(I2,M3)
            end
        end
    end

    for i ∈ 1:a

        if (i==1)
            if (i==c1)
                M4 = conket1
            elseif (i==c2)
                M4 = conket0
            else
                M4 = I2
            end
        else
            if (i==c1)
                M4 = kron(conket1,M4)
            elseif (i==c2)
                M4 = kron(conket0,M4)
            else
                M4 = kron(I2,M4)
            end
        end
    end
    return push!(C.gates,M1+M2+M3+M4)

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

function shots(C::Circuit, nshots::Int)
    a    
end