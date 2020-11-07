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

function CU!(C::Circuit,U::AbstractSparseArray{ComplexF64,Int},control::Int,target::Int)
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

function CU!(C::Circuit,λ::Float64,control::Int,target::Int)
    #Implentation of CU gate using Kronnecker products.
    a = length(C.qreg)
    U = sparse([1.0  0.0
                0.0  exp(im*λ)])

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
    elseif ((maximum(controls)>a) || (maximum(targets)>a)) || ((minimum(controls)<=0) || (minimum(targets)<=0))
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
    if ((control>a) || (maximum(targets)>a)) || ((control<=0) || (minimum(targets)<=0))
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
    if ((maximum(controls)>a) || (target>a)) || ((minimum(controls)<=0) || (target<=0))
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
    elseif ((maximum(controls)>a) || (maximum(targets)>a)) || ((minimum(controls)<=0) || (minimum(targets)<=0))
        error("The target and control qubits must all lie within the register!")
    end

    for i in 1:b
        CU!(C,U,controls[i],targets[i])
    end
end

function CU!(C::Circuit,U::AbstractSparseArray{ComplexF64,Int},control::Int,targets::AbstractVector{Int})
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    if ((control>a) || (maximum(targets)>a)) || ((control<=0) || (minimum(targets)<=0))
        error("The target and control qubits must all lie within the register!")
    end

    for i in targets
        CU!(C,U,control,i)
    end
end

function CU!(C::Circuit,U::AbstractSparseArray{ComplexF64,Int},controls::AbstractVector{Int},target::Int)
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    if ((maximum(controls)>a) || (target>a)) || ((minimum(controls)<=0) || (target<=0))
        error("The target and control qubits must all lie within the register!")
    end

    for i in controls
        CU!(C,U,i,target)
    end
end

function CU!(C::Circuit,λ::Float64,controls::AbstractVector{Int},targets::AbstractVector{Int})
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    b = length(controls)
    if length(targets)!==b
        error("The set of controls and targets must be the same length or the operation is ambiguous.")
    elseif ((maximum(controls)>a) || (maximum(targets)>a)) || ((minimum(controls)<=0) || (minimum(targets)<=0))
        error("The target and control qubits must all lie within the register!")
    end

    for i in 1:b
        CU!(C,λ,controls[i],targets[i])
    end
end

function CU!(C::Circuit,λ::Float64,control::Int,targets::AbstractVector{Int})
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    if ((control>a) || (maximum(targets)>a)) || ((control<=0) || (minimum(targets)<=0))
        error("The target and control qubits must all lie within the register!")
    end

    for i in targets
        CU!(C,λ,control,i)
    end
end

function CU!(C::Circuit,λ::Float64,controls::AbstractVector{Int},target::Int)
    #convenience method for doing mutliple CU operations at a time.
    a = length(C.qreg)
    if ((maximum(controls)>a) || (target>a)) || ((minimum(controls)<=0) || (target<=0))
        error("The target and control qubits must all lie within the register!")
    end

    for i in controls
        CU!(C,λ,i,target)
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