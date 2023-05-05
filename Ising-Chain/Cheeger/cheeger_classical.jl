using PythonCall
plt = pyimport("matplotlib.pyplot")

include("../ising_hamiltonian.jl")
include("../../classical_P.jl")
include("../../conductance.jl")





beta = 10
N = 4
couplings = ones(N)
couplings[end] = 0
h = 0
Hc = ising_ham(N, couplings, h)
energies = [Hc[i,i] for i in (1:2^N)]
println("Energies ", energies)
P_MH = MH_local(N,Hc, beta)

phi_min, min_A, all_phi_A, all_subset = min_conductance_search(Hc, P_MH, N,beta)
println("Min set ", min_A)
println("Min set energies ", [Hc[i,i] for i in min_A])
println("phi min ", phi_min)

# indicies = []
# for (i,phi) in pairs(all_phi_A)
#     if phi <(phi_min+0.01*phi_min)
#         append!(indicies,i)
#         println(all_subset[i], " with phi ", phi)
#     end
# end
println(conductance(Hc, P_MH, N, [1],beta))



# # Scaling with N 
beta = 5
N_values = (4:11)
phi_all_uniform = zeros(length(N_values))
phi_all_local = zeros(length(N_values))
for (i,N) in pairs(N_values)
    N = N_values[i]
    println("working on N ",N)
    couplings = ones(N)
    couplings[end] = 0
    Hc = ising_ham(N, couplings, 0)
    P_MH_uniform = MH_uniform(N,Hc, beta)
    phi_min = conductance(Hc, P_MH_uniform, N, [1],beta)
    phi_all_uniform[i] = phi_min
    P_MH_loc = MH_local(N,Hc, beta)
    # phi_min = conductance(Hc, P_MH_loc, N, (2:2^N-1),beta)
    phi_min = conductance(Hc, P_MH_loc, N, [1],beta)
    phi_all_local[i] = phi_min
end
save_object("Data/Ising-Chain/Classical/CheegerMHUniformScaling", phi_all_uniform)
save_object("Data/Ising-Chain/Classical/CheegerMHLocScaling", phi_all_local)



