using Statistics
using GLM
using LaTeXStrings
using PythonCall
plt = pyimport("matplotlib.pyplot")
include("ising_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")


function JW_energy(N, h, J1, J2)
    if N % 2 == 0
        ks_0 = (2*pi/N)*range(1,(N)/2) .- pi/N
        ks_1 = (2*pi/N)*range(1,(N)/2-1) 
    else
        println("Requires N to be even")
        return
    end
    Ek_0 = []
    for k in ks_0
        E = 2*sqrt((-h-J1*cos(k)-J2*cos(2*k))^2 +(J1*sin(k)+J2*sin(2*k))^2)
        push!(Ek_0, E)
    end
    Ek_1 = []
    for k in ks_1
        E = 2*sqrt((-h-J1*cos(k)-J2*cos(2*k))^2 +(J1*sin(k)+J2*sin(2*k))^2)
        push!(Ek_1, E)
    end
    return Ek_0, Ek_1
end


# Unit Test 
beta = 5
N = 8
couplings = ones(N)
couplings[end] = 0
H = ising_ham(N, couplings, 0)
energies = [H[i,i] for i in (1:2^N)]
# println("energies: ", energies)

P_glab = glab_local(N,H, beta)
eg,vg = eigen(Matrix(transpose(P_glab)))
gapg = abs(1-abs(eg[end-1]))
println("Gap glaub: ", gapg)
println("Gap second glaub: ", abs(1-abs(eg[end-2])))

h = 0.25*(1+1/cosh(2*beta))
J1 = 0.5*(tanh(2*beta))
J2 = 0.25*(1-1/cosh(2*beta))
Ek_0,Ek_1 = JW_energy(N, h, J1, J2)
println(Ek_0)
println(Ek_1)
E_0 = -sum(Ek_0)+N/2
E_1 = -sum(Ek_1)-2*J1+N/2
println("Jordan Wigner Gap: ", (E_1-E_0)/N)
println("Jordan Wigner Second Gap: ", (2*Ek_0[end])/N)




# N_values = (4:2:100)
# num = length(N_values)
# gap_glab_loc = zeros(num)
# gap2_glab_loc = zeros(num)
# beta = 5
# h = 0.25*(1+1/cosh(2*beta))
# J1 = 0.5*(tanh(2*beta))
# J2 = 0.25*(1-1/cosh(2*beta))
# for j in (1:num)
#     N = N_values[j]
#     Ek_0,Ek_1 = JW_energy(N, h, J1, J2)
#     E_0 = -sum(Ek_0)+N/2
#     E_1 = -sum(Ek_1)-2*J1+N/2
#     gap_glab_loc[j] = (E_1-E_0)/N
#     gap2_glab_loc[j]= (2*Ek_0[end])/N
# end
# save_object("Data/Ising-Chain/Classical/GlaubLocScalingAnalytic",gap2_glab_loc)

# x = log.(N_values)
# y = log.(gap_glab_loc) 
# fit  = lm(@formula(y ~ x), (;x,y))
# display(fit)
# param = coef(fit)
# scale = -round(param[2], digits=2)

# y = log.(gap2_glab_loc) 
# fit2  = lm(@formula(y ~ x), (;x,y))
# display(fit2)
# param2 = coef(fit2)
# scale2 = -round(param2[2], digits=2)

# label_scatter = L"$\delta_0$ Glauber Local Proposal $N^{-b}$ b= "*string(scale)
# plt.scatter(N_values,  gap_glab_loc, color = "tab:blue")
# x = range(4,50, length= 1000)
# plt.plot(x,[exp(param[2]*log(i)+param[1]) for i in x],linestyle = "dashed", label = label_scatter,color = "tab:blue")
# plt.ylabel(L"$\delta$")

# label_scatter2 = L"$\delta_1$ Glauber Local Proposal $N^{-b}$ b= "*string(scale2)
# plt.scatter(N_values,  gap2_glab_loc, color = "tab:red")
# plt.plot(x,[exp(param2[2]*log(i)+param2[1]) for i in x],linestyle = "dashed", label = label_scatter2,color = "tab:red")
# plt.ylabel(L"$\delta$")

# plt.xlabel(L"$N$")
# plt.yscale("log")
# plt.xscale("log")
# plt.grid(which="major", zorder=-1.0)
# plt.legend()
# name = "Figures/Ising-Chain/glauberAnalyticBeta"*string(beta)*".svg"
# plt.savefig(name)
# plt.show()
