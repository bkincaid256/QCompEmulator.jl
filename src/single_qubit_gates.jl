######Single qubit gates ###########

function U!(C::Circuit, U::AbstractArray{ComplexF64,2}, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    if (bitpos>a) || (bitpos<=0)
        error("The bit specified must be within the register of the circuit!")
    end
    Uf = sparse(U)
    M = similar(Uf)
    I2 = sparse(I*(1.0+0.0*im),2,2)
    for  i ∈ 1:a
        if (i == 1) && (i !== bitpos)
            copy!(M,I2)
        elseif (i == 1) && (i == bitpos)
            copy!(M,Uf)
        else  
            if (i == bitpos)
                M = kron(Uf,M)
            else
                M = kron(I2,M)
            end
        end
    end
   return push!(C.gates, M)
end

function U!(C::Circuit, U::AbstractSparseArray{ComplexF64,Int}, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    a = length(C.qreg)
    if (bitpos>a) || (bitpos<=0)
        error("The bit specified must be within the register of the circuit!")
    end
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

function U!(C::Circuit, λ::Float64, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#

    U = sparse([1.0  0.0
                0.0  exp(im*λ)])
    a = length(C.qreg)
    if (bitpos>a) || (bitpos<=0)
        error("The bit specified must be within the register of the circuit!")
    end
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

##Specialized methods for doing multiple gates sequentially
function U!(C::Circuit, U::AbstractArray{ComplexF64,2}, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        U!(C,U,i)
    end
end

function U!(C::Circuit, U::AbstractSparseArray{ComplexF64,Int}, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        U!(C,U,i)
    end
end

##Specialized methods for doing multiple gates sequentially
function U!(C::Circuit, λ::Float64, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        U!(C,λ,i)
    end
end

function Udagger!(C::Circuit, λ::Float64, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#

    U = sparse([1.0  0.0
                0.0  conj(exp(im*λ))])
    a = length(C.qreg)
    if (bitpos>a) || (bitpos<=0)
        error("The bit specified must be within the register of the circuit!")
    end
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

##Specialized methods for doing multiple gates sequentially
function Udagger!(C::Circuit, λ::Float64, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    for i in set
        Udagger!(C,λ,i)
    end
end
function H!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    H = sparse([1.0/sqrt(2)+0.0*im  1.0/sqrt(2)+0.0*im
                1.0/sqrt(2)+0.0*im -1.0/sqrt(2)+0.0*im])
    U!(C,H,bitpos)
end

function H!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    H = sparse([1.0/sqrt(2)+0.0*im  1.0/sqrt(2)+0.0*im
                1.0/sqrt(2)+0.0*im -1.0/sqrt(2)+0.0*im])
    U!(C,H,set)
end

function T!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    T = sparse([1.0+0.0*im  0.0+0.0*im
                0.0+0.0*im  exp(im*π/4)]) 
    U!(C,T,bitpos)
end

function T!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    T = sparse([1.0+0.0*im  0.0+0.0*im
                0.0+0.0*im  exp(im*π/4)]) 
    U!(C,T,set)
end

######## Pauli Gates ##########
function X!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    X = sparse([0.0+0.0*im  1.0+0.0*im
                1.0+0.0*im  0.0+0.0*im]) 
    U!(C,X,bitpos)
end

function X!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    X = sparse([0.0+0.0*im  1.0+0.0*im
                1.0+0.0*im  0.0+0.0*im]) 
    U!(C,X,set)
end

function Y!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    Y = sparse([0.0+0.0*im  0.0-1.0*im
                0.0+1.0*im  0.0+0.0*im]) 
    U!(C,Y,bitpos)
end

function Y!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    Y = sparse([0.0+0.0*im  0.0-1.0*im
                0.0+1.0*im  0.0+0.0*im]) 
    U!(C,Y,set)
end

function Z!(C::Circuit, bitpos::Int)
    #=Implementation that adds the new matrix into the gates field 
    of the circuit=#
    Z = sparse([1.0+0.0*im  0.0+0.0*im
                0.0+0.0*im -1.0+0.0*im]) 
    U!(C,Z,bitpos)
end

function Z!(C::Circuit, set::AbstractVector{Int})
    #=Implementation that gives user ability to define 
    a range of qubits to affect with the specified gate=#
    Z = sparse([1.0+0.0*im  0.0+0.0*im
                0.0+0.0*im -1.0+0.0*im]) 
    U!(C,Z,set)
end