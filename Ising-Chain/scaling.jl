using Arpack
using JLD2

include("ising_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")

h = 0
N_values = (4:11) 
num = length(N_values)
beta = 5

# Get best gap 
# num_values = 300
# eta_best = zeros(length(N_values))
# alpha_best = zeros(length(N_values))
# for i in (1:length(N_values))
#     N = N_values[i]
#     name = "Data/Ising/qHMC/GridSearch/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
#     gap_all= load_object(name)
#     max_param = 50
#     max,cord = findmax(gap_all)
#     eta_best[i] = (cord[2]-1)/num_values*max_param
#     alpha_best[i] = (cord[1]-1)/num_values*max_param
# end


gap_glab = zeros(num)
gap_glab_loc = zeros(num)
gap_MH = zeros(num)
gap_MH_loc = zeros(num)
# gap_rand = zeros(num)
# gap_qHMC = zeros(num)
for j in (1:num)
    N = N_values[j]
    couplings = ones(N)
    # couplings[end] = 0
    println(" Working on N = ",N)
    H = ising_ham(N, couplings, h)

    P_glab_uniform = glab_uniform(N,H, beta)
    e,v  = eigs(P_glab_uniform, nev = 2, which=:LM)
    gap_glab[j] = abs(1-abs(e[2]))

    P_glab_loc = glab_local(N,H, beta)
    e,v  = eigs(P_glab_loc, nev = 2, which=:LM)
    gap_glab_loc[j] = abs(1-abs(e[2]))

    P_MH_uniform = MH_uniform(N,H, beta)
    e,v  = eigs(P_MH_uniform, nev = 2, which=:LM)
    gap_MH[j] = abs(1-abs(e[2]))

    P_MH_loc = MH_local(N,H, beta)
    e,v  = eigs(P_MH_loc, nev = 2, which=:LM)
    gap_MH_loc[j] = abs(1-abs(e[2]))

    # Hm = mixing_ham(N) 
    # M_qHMC = mixing_qHMC(N, H,alpha_best[j],eta_best[j],beta, Hm)
    # try
    #     e,v  = eigs(M_qHMC, nev = 2, which=:LM)
    #     gap_qHMC[j] = abs(1-abs(e[2]))
    # catch
    #     e,v = eigen(Matrix(M_qHMC ))
    #     gap_qHMC[j]= abs(1-abs(e[end-1]))
    # end
    # M_rand = random_mixing(N,H,beta,Hm, 100)
    # try
    #     e,v  = eigs(M_rand, nev = 2, which=:LM)
    #     gap_rand[j] = abs(1-abs(e[2]))
    # catch
    #     e,v = eigen(Matrix(M_rand))
    #     gap_rand[j]= abs(1-abs(e[end-1]))
    # end
end
save_object("Data/Ising-Chain/Classical/GlaubScaling", gap_glab)
save_object("Data/Ising-Chain/Classical/GlaubLocScaling", gap_glab_loc)
# save_object("Data/Ising-Chain/Classical/MHUniformScaling", gap_MH)
# save_object("Data/Ising-Chain/Classical/MHLocScaling", gap_MH_loc)
# save_object("Data/Ising/qHMC/RandScaling", gap_rand)
# save_object("Data/Ising/qHMC/BestScaling", gap_qHMC)
