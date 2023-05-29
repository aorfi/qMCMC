using Arpack
using JLD2
using Random, Distributions
using Statistics
using StatsBase
Random.seed!(123)

include("SK_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")


N_values = (4:10) 
num = length(N_values)
beta = 5
runs = 100


# gap_all_MH = zeros(num,2)
# gap_all_MH_loc = zeros(num,2)

# for i in (1:num)
#     N = N_values[i]
#     println(" Working on N = ",N)
#     gap_MH = zeros(runs)
#     gap_MH_loc = zeros(runs)
#     for j in (1:runs)
#         println(" Working on run = ",j)
#         h = rand(Normal(0,1),N)
#         couplings = rand(Normal(0,1),sum(1:N-1))
#         H = SK_ham(N,couplings,h)
#         P_MH_uniform = MH_uniform(N,H, beta)
#         try
#             e,v  = eigs(P_MH_uniform, nev = 2, which=:LM)
#             gap_MH[j] = abs(1-abs(e[2]))
#         catch
#             e,v = eigen(Matrix(P_MH_uniform))
#             gap_MH[j]= abs(1-abs(e[end-1]))
#         end
#         P_MH_loc = MH_local(N,H, beta)
#         try
#             e,v  = eigs(P_MH_loc, nev = 2, which=:LM)
#             gap_MH_loc[j] = abs(1-abs(e[2]))
#         catch
#             e,v = eigen(Matrix(P_MH_loc))
#             gap_MH_loc[j]= abs(1-abs(e[end-1]))
#         end
#     end
#     gap_all_MH[i,1] = mean(gap_MH)
#     gap_all_MH[i,2] = stdm(gap_MH, mean(gap_MH))
#     gap_all_MH_loc[i,1] = mean(gap_MH_loc)
#     gap_all_MH_loc[i,2] = stdm(gap_MH_loc, mean(gap_MH_loc))
# end
# save_object("Data/SK/Classical/MHScaling", gap_all_MH)
# save_object("Data/SK/Classical/MHLocScaling", gap_all_MH_loc)




num_q = length(N_values[1:end-2])
gap_av = zeros(num_q,2)
gap_best = zeros(num_q,2)

num_values = 50
kappa_values = range(0,50, length=num_values)
eta_values = range(0,50, length=num_values)


for i in (1:num_q)
    N = N_values[i]
    name = "Data/SK/qMCMC/Grid-Search/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
    gap_all= load_object(name)
    gap_best_runs = zeros(runs)
    gap_av_runs = zeros(runs)
    for j in (1:runs)   
        log_data = -log.(1 .- gap_all[:,:,j])
        max,cord = findmax(gap_all[:,:,j])
        gap_best_runs[j] = max
        gap_av_runs[j] = mean(log_data)
    end
    gap_av[i,1] = mean(gap_av_runs)
    gap_av[i,2] = stdm(gap_av_runs, mean(gap_av_runs))
    # println(stdm(gap_av_runs, mean(gap_av_runs)))
    gap_best[i,1] = mean(gap_best_runs)
    gap_best[i,2] = stdm(gap_best_runs, mean(gap_best_runs))
end
save_object("Data/SK/qMCMC/AverageGapScaling", gap_av)
save_object("Data/SK/qMCMC/LargestGapScaling", gap_qHMC)
println("saved")


