using Arpack
using JLD2
using Random, Distributions
Random.seed!(123)

include("SK_hamiltonian.jl")
include("../quantum_P.jl")


runs = 2
N_values = (4:4) 
for i in (1:length(N_values))
    N = N_values[i]
    beta = 5
    num_values = 10
    kappa_values = range(0,50, length=num_values)
    eta_values = range(0,50, length=num_values)
    gap_all = zeros(num_values,num_values,runs)
    for kappa_i in (1:num_values)
        println("kappa: ", kappa_i)
        for eta_i in (1:num_values)
            println("  eta: ", eta_i)
            Hm = mixing_ham(N) 
            for j in (1:runs)
                h = rand(Normal(0,1),N)
                couplings = rand(Normal(0,1),sum(1:N-1))
                H = SK_ham(N,couplings,h)
                P = mixing_qHMC(N, H,kappa_values[kappa_i],eta_values[eta_i],beta, Hm) 
                e,v  = eigs(P, nev = 2, which=:LM)
                try
                    e,v  = eigs(P, nev = 2, which=:LM)
                    gap_all[kappa_i,eta_i,j]  = abs(1-abs(e[2]))
                catch
                    e,v = eigen(Matrix(P))
                    gap_all[kappa_i,eta_i,j] = abs(1-abs(e[end-1]))
                end
            end
        end
    end
    name = "Data/SK/qMCMC/Grid-Search/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
    save_object(name, gap_all)
end
