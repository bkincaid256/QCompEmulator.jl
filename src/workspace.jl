using  LinearAlgebra, StaticArrays, BenchmarkTools

mutable struct Qubit #type for holding qubit information
    bitpos::Int
    vec::Vector{ComplexF64}
end

mutable struct Circuit
    qreg::Array{Qubit,1}
    gates::Array{Matrix}
end

function makeregister(nqubits::Int)
    #function for generating a quantum register of qubits
    zerovec = @SVector [1.0+0.0*im, 0.0+0.0*im]
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

function ControlX!(C::Circuit,target::Integer,control::Integer)
    
    
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
        if (target-control) == 1    

            for i ∈ 0:a-2
                if (i == 0) && (control-1 == 0)
                    M = kron(I(1),cX)
                elseif (i == 0) && (control-1 !== 0)
                    M = kron(I(1),I(2))
                else
                    if (i == control-1)
                        M = kron(cX,M)
                    else
                        M = kron(I(2),M)
                    end
                end
            end
        elseif (control-target) == 1    

            for i ∈ 0:a-2
                if (i == 0) && (target-1 == 0)
                    M = kron(I(1),cX)
                elseif (i == 0) && (target-1 !== 0)
                    M = kron(I(1),I(2))
                else
                    if (i == target-1)
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