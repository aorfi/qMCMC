include("../ising_hamiltonian.jl")
include("../../classical_P.jl")
include("../../conductance.jl")
include("../../quantum_P.jl")


beta = 5
num_values = 10

N = 4

couplings = ones(N)
Hc = ising_ham(N, couplings, 0)
energies = [Hc[i,i] for i in (1:2^N)]
# println("Energies ", [Hc[i,i] for i in (1:2^N)])
s1 = [i for i in (1:2^N) if energies[i]== (energies[1]+4)]

all_subsets = collect(subsets((1:2^N)))[2:end]

kappa_values = range(0,50, length=num_values)
eta_values = range(0,50, length=num_values)

phi_all = zeros(num_values,num_values)
indexes = zeros(num_values,num_values)
for kappa_i in (1:num_values)
    println("kappa: ", kappa_i)
    for eta_i in (1:num_values)
        println("  eta: ", eta_i)
        Hm = mixing_ham(N) 
        P = mixing_qHMC(N, Hc,kappa_values[kappa_i],eta_values[eta_i],beta, Hm) 
        phi_min, min_A, all_phi_A,subset_small= min_conductance_search(Hc, P, N, beta)
        phi_all[kappa_i,eta_i]  = phi_min
    end
end
name = "Data/Ising-Chain/qMCMC/CheegerGrid"*string(num_values)*"N"*string(N)*"beta"*string(beta)
save_object(name, phi_all)