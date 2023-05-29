using Arpack
using JLD2
using Random, Distributions
using Statistics
Random.seed!(123)

include("SK_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")


N = 8
runs = 100
num_points = 50
temp = 10 .^ (range(-3,stop=3,length=num_points))
beta_values = 1 ./ temp

# Get classical data
gap_all_MH = zeros(num_points,2)
gap_all_MH_loc = zeros(num_points,2)

for (i,beta) in pairs(beta_values)
    println(" Working on T = ", temp[i])
    gap_MH = zeros(runs)
    gap_MH_loc = zeros(runs)
    for j in (1:runs)
        println(" Working on run = ",j)
        h = rand(Normal(0,1),N)
        couplings = rand(Normal(0,1),sum(1:N-1))
        H = SK_ham(N,couplings,h)
        # P_MH_uniform = MH_uniform(N,H, beta)
        # try
        #     e,v  = eigs(P_MH_uniform, nev = 2, which=:LM)
        #     gap_MH[j] = abs(1-abs(e[2]))
        # catch
        #     e,v = eigen(Matrix(P_MH_uniform))
        #     gap_MH[j]= abs(1-abs(e[end-1]))
        # end
        P_MH_loc = MH_local(N,H, beta)
        try
            e,v  = eigs(P_MH_loc, nev = 2, which=:LM)
            gap_MH_loc[j] = abs(1-abs(e[2]))
        catch
            e,v = eigen(Matrix(P_MH_loc))
            gap_MH_loc[j]= abs(1-abs(e[end-1]))
        end
    end
    # gap_all_MH[i,1] = mean(gap_MH)
    # gap_all_MH[i,2] = stdm(gap_MH, mean(gap_MH))
    gap_all_MH_loc[i,1] = mean(gap_MH_loc)
    gap_all_MH_loc[i,2] = stdm(gap_MH_loc, mean(gap_MH_loc))
end
# save_object("Data/SK/Classical/N"*string(N)*"MHTempAv", gap_all_MH)
save_object("Data/SK/Classical/N"*string(N)*"MHLocTempAv", gap_all_MH_loc)


## get quantum data
# for (i,beta) in pairs(beta_values)
#     println("Working on beta ", beta)
#     num_values = 2
#     kappa_values = range(0,50, length=num_values)
#     eta_values = range(0,50, length=num_values)
#     gap_all = zeros(num_values,num_values,runs)
#     for kappa_i in (1:num_values)
#         println("kappa: ", kappa_i)
#         for eta_i in (1:num_values)
#             println("  eta: ", eta_i)
#             Hm = mixing_ham(N) 
#             for j in (1:runs)
#                 println(" Working on run = ",j)
#                 h = rand(Normal(0,1),N)
#                 couplings = rand(Normal(0,1),sum(1:N-1))
#                 H = SK_ham(N,couplings,h)
#                 P = mixing_qHMC(N, H,kappa_values[kappa_i],eta_values[eta_i],beta, Hm) 
#                 e,v  = eigs(P, nev = 2, which=:LM)
#                 try
#                     e,v  = eigs(P, nev = 2, which=:LM)
#                     gap_all[kappa_i,eta_i,j]  = abs(1-abs(e[2]))
#                 catch
#                     e,v = eigen(Matrix(P))
#                     gap_all[kappa_i,eta_i,j] = abs(1-abs(e[end-1]))
#                 end
#             end
#         end
#     end
#     name = "Data/SK/qMCMC/Grid-Search/Temperature/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
#     save_object(name, gap_all)
# end

# # # # Average temperature gap
# gap_all_av = zeros(num_points,2)
# gap_all_best = zeros(num_points,2)
# num_values = 10
# for (i,beta) in pairs(beta_values)

#     name = "Data/SK/qMCMC/Grid-Search/Temperature/Small/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
#     gap_all= load_object(name)
#     gap_best = zeros(runs)
#     for j in (1:runs)
#         gap = gap_all[:,:,j]
#         max,cord = findmax(gap)
#         gap_best[j] = max
#     end
#     gap_all_av[i,1] = mean(gap_all)
#     gap_all_av[i,2] = std(gap_all)
#     gap_all_best[i,1] = mean(gap_best)
#     gap_all_best[i,2] = stdm(gap_best, mean(gap_best))
# end

# save_object("Data/SK/qMCMC/N"*string(N)*"AverageTemp", gap_all_av)
# save_object("Data/SK/qMCMC/N"*string(N)*"BestTemp", gap_all_best)

