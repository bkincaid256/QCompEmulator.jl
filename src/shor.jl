using LaTeXStrings
include("workspace.jl")
pyplot()

N = 27
a = 4

#Calculate the plotting data 
xvals = 0:N-1
yvals = [mod(a^x,N) for x ∈ xvals]

plot(xvals, yvals, linewidth=1, ls = :dot, marker = :x,
xlabel = "\$x\$", ylabel = "\$$a^x\$ mod \$$N\$",
title = "Periodicity of modular exponentiation", legend = false, grid = false, mc = :red, lc = :red)
xticks!(0:5:xvals[end]+1)
gui()

r = 0
i = 2
while i < length(xvals)
    if yvals[1] == yvals[i]
        return r = i-1
        break
    end
    if i == length(xvals)
        println("No period was found")
        break
    end
    i += 1
end
plot!([0:1:r,r:-1:0],ones(length(0:1:r),2),
arrow = :arrow, lc = :black)
annotate!((r/2,1.5, Plots.text("\$r = $r\$", 8,:black,:center)))

##
function controlbits(a::Int, n::Int)
    # a is the starting random integer, n is the number of qubits in the circuit
    digits(a,base=2,pad=n) 
end

function phaseadd!(C::Circuit, controlbits::AbstractArray{Int})
    n = length(controlbits)
    for m in 1:n
        for i in m:n
            if controlbits[n-i+1] == 1
                k = i-m+1
                λ = 2π/(2^(k))
                U!(C,λ,n-m+1)
            end
        end
    end

   # for qubit in 1:div(n,2)-1
   #     Swap!(C, qubit, n-qubit)
   # end
end

function phaseadd_dagger!(C::Circuit, controlbits::AbstractArray{Int})
    n = length(controlbits)
    for m in n:-1:1
        for i in n:-1:m
            if controlbits[n-i+1] == 1
                λ = 2π/2^(i-m+1)
                U!(C,-λ,n-m+1)
            end
        end
    end
end

#for m in 1:n
#    for i in m:n
#        println(n-i+1," ",i-m+1)
#    end
#end
#
#for m in 5:-1:1
#    for i in 5:-1:m
#        println(5-i+1," ",i-m+1)
#    end
#end

function qft_dagger!(C::Circuit, n::Int)
    #Don't forget the swaps!
    #for qubit in 1:div(n,2)-1
    #    Swap!(C, qubit, n-qubit)
    #end
    for j in 1:n-1
        for m in 1:j-1
            CU!(C,-π/(2^(j-m)), m, j)
        end
        H!(C,j)
    end
end

function qft!(C::Circuit, n::Int)
    #Don't forget the swaps!
    for j in n-1:-1:1
        H!(C,j)
        for m in j-1:-1:1
            CU!(C,π/(2^(j-m)), m, j)
        end
    end
    #for qubit in 1:div(n,2)-1
    #    Swap!(C, qubit, n-qubit)
    #end
end

n = 3
c1 = Circuit(n)
a = 3

#X!(c1,1)
#X!(c1,2)
controls = controlbits(a,n-1)
qft!(c1,n)
phaseadd!(c1,controls)
#phaseadd_dagger!(c1,controls)
qft_dagger!(c1,n)
#qft_dagger!(c1,n)
#qft!(c1,n)
dat = shots(c1, 10_000)
graph!(dat)