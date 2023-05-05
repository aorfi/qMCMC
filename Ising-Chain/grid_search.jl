using Arpack
using JLD2

include("ising_hamiltonian.jl")
include("../quantum_P.jl")

N_values = (4:4) 
for i in (1:length(N_values))
    N = N_values[i]
    h = 0
    couplings = ones(N)
    # couplings[end] = 0
    H = ising_ham(N, couplings, h)
    beta = 5
    num_values = 300
    kappa_values = range(0,50, length=num_values)
    eta_values = range(0,50, length=num_values)

    gap_all = zeros(num_values,num_values)
    for kappa_i in (1:num_values)
        println("kappa: ", kappa_i)
        for eta_i in (1:num_values)
            println("  eta: ", eta_i)
            Hm = mixing_ham(N) 
            P = mixing_qHMC(N, H,kappa_values[kappa_i],eta_values[eta_i],beta, Hm) 
            e,v  = eigs(P, nev = 2, which=:LM)
            gap_all[kappa_i,eta_i]  = abs(1-abs(e[2]))   
        end
    end
    name = "Data/Ising-Chain/qMCMC/Grid-Search/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
    save_object(name, gap_all)
end
