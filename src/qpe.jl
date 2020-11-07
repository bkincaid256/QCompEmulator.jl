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
    end
    reps *= 2
end

function qft_dagger!(C::Circuit, n::Int)
    #Don't forget the swaps!
    for qubit in 1:div(n,2)-1
        Swap!(C, qubit, n-qubit)
    end
    for j in 1:n-1
        for m in 1:j-1
            CU!(C,-π/(2^(j-m)), m, j)
        end
        H!(C,j)
    end
end

qft_dagger!(qpe,4)
data = shots(qpe,10_000)
graph!(data)

#qpe2
qpe2 = Circuit(4)

X!(qpe2,4)
H!(qpe2,1:3)

λ2 = 2*π/3
reps = 1
for counting_qubit ∈ 1:3
    for i ∈ 1:reps
        CU!(qpe2, λ2, counting_qubit, 4)
    end
    reps *= 2
end

qft_dagger!(qpe2,4)
data2 = shots(qpe2,10_000)
graph!(data2)