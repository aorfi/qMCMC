using Arpack
using JLD2
using Statistics

include("ising_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")

h = 0
N_values = (4:14) 
num = length(N_values)
beta = 5


# gap_glab = zeros(num)
# gap_glab_loc = zeros(num)
# gap_MH = zeros(num)
# gap_MH_loc = zeros(num)
# for j in (1:num)
#     N = N_values[j]
#     couplings = ones(N)
#     # couplings[end] = 0
#     println(" Working on N = ",N)
#     H = ising_ham(N, couplings, h)

#     # P_glab_uniform = glab_uniform(N,H, beta)
#     # e,v  = eigs(P_glab_uniform, nev = 2, which=:LM)
#     # gap_glab[j] = abs(1-abs(e[2]))

#     # P_glab_loc = glab_local(N,H, beta)
#     # e,v  = eigs(P_glab_loc, nev = 2, which=:LM)
#     # gap_glab_loc[j] = abs(1-abs(e[2]))

#     P_MH_uniform = MH_uniform(N,H, beta)
#     e,v  = eigs(P_MH_uniform, nev = 2, which=:LM)
#     gap_MH[j] = abs(1-abs(e[2]))

#     P_MH_loc = MH_local(N,H, beta)
#     e,v  = eigs(P_MH_loc, nev = 2, which=:LM)
#     gap_MH_loc[j] = abs(1-abs(e[2]))

# end
# # save_object("Data/Ising-Chain/Classical/GlaubScaling", gap_glab)
# # save_object("Data/Ising-Chain/Classical/GlaubLocScaling", gap_glab_loc)
# save_object("Data/Ising-Chain/Classical/MHUniformScaling", gap_MH)
# save_object("Data/Ising-Chain/Classical/MHLocScaling", gap_MH_loc)


num_q = length(N_values[1:end-6])
gap_av = zeros(num_q,2)
gap_qHMC = zeros(num_q)

num_values = 300
kappa_values = range(0,50, length=num_values)
eta_values = range(0,50, length=num_values)


for i in (1:num_q )
    N = N_values[i]
    name = "Data/Ising-Chain/qMCMC/Grid-Search/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
    gap_all= load_object(name)
    max,cord = findmax(gap_all)
    gap_qHMC[i] = max
    gap_av[i,1] = mean(gap_all)
    gap_av[i,2] = stdm(gap_all, mean(gap_all))
end
save_object("Data/Ising-Chain/qMCMC/AverageGapScaling", gap_av)
save_object("Data/Ising-Chain/qMCMC/LargestGapScaling", gap_qHMC)