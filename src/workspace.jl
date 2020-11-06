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

############# Controlled Gates and 2 qubit Gates#############
function CU!(C::Circuit, U::AbstractArray{ComplexF64,2},control::Int,target::Int)
    #Implentation of CU gate using Kronnecker products.
    Uf = sparse(U)

    a = length(C.qreg)

    if (target==control)
        error("Target qubit and control qubit must be Different!")
    elseif ((target>a) || (control>a)) || ((target<=0) || (control<=0)) 
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
                M1 = Uf
            else
                M1 = I2
            end
        else
            if (i==control)
                M1 = kron(conket1,M1)
            elseif (i==target)
                M1 = kron(Uf,M1)
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

function CU!(C::Circuit, U::AbstractSparseArray{ComplexF64,Int},control::Int,target::Int)
    #Implentation of CU gate using Kronnecker products.
    a = length(C.qreg)

    if (target==control)
        error("Target qubit and control qubit must be Different!")
    elseif ((target>a) || (control>a)) || ((target<=0) || (control<=0)) 
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
                M1 = U
            else
                M1 = I2
            end
        else
            if (i==control)
                M1 = kron(conket1,M1)
            elseif (i==target)
                M1 = kron(U,M1)
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
###Methods for doing a number of gates in a row
function CU!(C::Circuit,U::AbstractArray{ComplexF64,2},controls::AbstractVector{Int},targets::AbstractVector{Int})
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    b = length(controls)
    if length(targets)!==b
        error("The set of controls and targets must be the same length or the operation is ambiguous.")
    elseif ((max(controls)>a) || (max(targets)>a)) || ((min(controls)<=0) || (min(targets)<=0))
        error("The target and control qubits must all lie within the register!")
    end

    Uf = sparse(U)
    for i in 1:b
        CU!(C,Uf,controls[i],targets[i])
    end
end

function CU!(C::Circuit,U::AbstractArray{ComplexF64,2},control::Int,targets::AbstractVector{Int})
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    if ((control>a) || (max(targets)>a)) || ((control<=0) || (min(targets)<=0))
        error("The target and control qubits must all lie within the register!")
    end

    Uf = sparse(U)
    for i in targets
        CU!(C,Uf,control,i)
    end
end

function CU!(C::Circuit,U::AbstractArray{ComplexF64,2},controls::AbstractVector{Int},target::Int)
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    if ((max(controls)>a) || (target>a)) || ((min(controls)<=0) || (target<=0))
        error("The target and control qubits must all lie within the register!")
    end

    Uf = sparse(U)
    for i in controls
        CU!(C,Uf,i,target)
    end
end

function CU!(C::Circuit,U::AbstractSparseArray{ComplexF64,Int},controls::AbstractVector{Int},targets::AbstractVector{Int})
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    b = length(controls)
    if length(targets)!==b
        error("The set of controls and targets must be the same length or the operation is ambiguous.")
    elseif ((max(controls)>a) || (max(targets)>a)) || ((min(controls)<=0) || (min(targets)<=0))
        error("The target and control qubits must all lie within the register!")
    end

    for i in 1:b
        CU!(C,U,controls[i],targets[i])
    end
end

function CU!(C::Circuit,U::AbstractSparseArray{ComplexF64,Int},control::Int,targets::AbstractVector{Int})
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    if ((control>a) || (max(targets)>a)) || ((control<=0) || (min(targets)<=0))
        error("The target and control qubits must all lie within the register!")
    end

    for i in targets
        CU!(C,U,control,i)
    end
end

function CU!(C::Circuit,U::AbstractSparseArray{ComplexF64,Int},controls::AbstractVector{Int},target::Int)
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    if ((max(controls)>a) || (target>a)) || ((min(controls)<=0) || (target<=0))
        error("The target and control qubits must all lie within the register!")
    end

    for i in controls
        CU!(C,U,i,target)
    end
end

function CX!(C::Circuit,control::Int,target::Int)
    #Implentation of CNOT gate using Kronnecker products.
    X = sparse([0.0+0.0*im 1.0+0.0*im;
                1.0+0.0*im 0.0+0.0*im])

    CU!(C,X,control,target)
end

function CX!(C::Circuit,controls::AbstractVector{Int},targets::AbstractVector{Int})
    #convenience method for doing mutliple CX operations at a time.
    X = sparse([0.0+0.0*im 1.0+0.0*im;
                1.0+0.0*im 0.0+0.0*im])

    CU!(C,X,controls,targets)
end

function CX!(C::Circuit,control::Int,targets::AbstractVector{Int})
    #convenience method for doing mutliple CX operations at a time.
    X = sparse([0.0+0.0*im 1.0+0.0*im;
                1.0+0.0*im 0.0+0.0*im])

    CU!(C,X,control,targets)
end

function CX!(C::Circuit,controls::AbstractVector{Int},target::Int)
    #convenience method for doing mutliple CX operations at a time.
    X = sparse([0.0+0.0*im 1.0+0.0*im;
                1.0+0.0*im 0.0+0.0*im])

    CU!(C,X,controls,target)
end

function CZ!(C::Circuit,control::Int,target::Int)
    #Implentation of CZ gate using Kronnecker products.
    Z = sparse([1.0+0.0*im  0.0+0.0*im;
                0.0+0.0*im -1.0+0.0*im])

    CU!(C,Z,control,target)
end

function CZ!(C::Circuit,controls::AbstractVector{Int},targets::AbstractVector{Int})
    #convenience method for doing mutliple CZ operations at a time.
    Z = sparse([1.0+0.0*im  0.0+0.0*im;
                0.0+0.0*im -1.0+0.0*im])

    CU!(C,Z,controls,targets)
end

function CZ!(C::Circuit,control::Int,targets::AbstractVector{Int})
    #convenience method for doing mutliple CZ operations at a time.
    Z = sparse([1.0+0.0*im  0.0+0.0*im;
                0.0+0.0*im -1.0+0.0*im])

    CU!(C,Z,control,targets)
end

function CZ!(C::Circuit,controls::AbstractVector{Int},target::Int)
    #convenience method for doing mutliple CZ operations at a time.
    Z = sparse([1.0+0.0*im  0.0+0.0*im;
                0.0+0.0*im -1.0+0.0*im])

    CU!(C,Z,controls,target)
end

function Swap!(C::Circuit,b1::Int,b2::Int)
    #quick and dirty implmentation of Swap Gate
    a = length(C.qreg)
    if b1==b2
        error("This ins't really a swap, you just tried to swap a qubit with itself!")
    elseif ((b1>a) || (b2>a)) || ((b1<=0) || (b2<=0))
        error("Both bits must be within the size of the register!")
    end
    CX!(C,b2,b1)
    CX!(C,b1,b2)
    CX!(C,b2,b1)
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
    
    finalket = Run_circ(C)
    ntrials = rand(nshots)
    ensemble = zeros(Int,length(finalket))
    trials!(ensemble, finalket, ntrials)
end

function trials!(ensemble::AbstractVector{Int}, ket::AbstractVector{ComplexF64}, probs::AbstractVector{Float64})
    nshots = length(probs)

    @inbounds for j ∈ 1:nshots
        local p = probs[j]
        local sumold = 0.0
        local sumnew = 0.0 
        
        @inbounds for i ∈ 1:length(ket)
            sumnew = sumold + real(conj(ket[i])*ket[i])
            if (p<=sumnew)
                ensemble[i] += 1
                break
            end
            sumold = sumnew
        end
    end
    return ensemble
end

function graph(data::AbstractVector)
    pyplot()
    max = length(data)-1
    nqubits = convert(Int,log(2,max+1))
    xs = [string(i,base=2,pad=nqubits) for i ∈ 0:max]

    bar(xs,data, leg=false, title="Simulation of Circuit", xlabel="qubit readout", ylabel="Counts", dpi = 300)
end