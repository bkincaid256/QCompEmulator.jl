include("workspace.jl")

#make circuit
qpe = Circuit(4)
X!(qpe,4)
H!(qpe,1:3)
λ = π/4

reps = 1
for counting_qubit ∈ 1:3
    for i ∈ 1:reps
        CU!(qpe, λ, counting_qubit, 4)
        reps *= 2
    end
end

function qft_dagger(C::Circuit, n::Int)
    #Don't forget the swaps!
   # for qubit in 1:div(n,2)-1
   #     Swap!(C, qubit, n-qubit-1)
   # end
    for j in 1:n-1
        for m in 1:j-1
            CU!(C,-π/(2^(j-m)), m, j)
        end
        H!(C,j)
    end
end
Swap!(qpe,1,3)

qft_dagger(qpe, 4)
