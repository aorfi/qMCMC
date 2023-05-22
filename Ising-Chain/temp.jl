using JLD2
using Statistics

include("ising_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")

N = 8
h = 0
couplings = ones(N)
# couplings[end] = 0
H = ising_ham(N, couplings, h)

num_points = 50
temp = 10 .^ (range(-3,stop=3,length=num_points))
beta_values = 1 ./ temp
gap_MH = zeros(num_points)
gap_MH_loc = zeros(num_points)
for (j,beta) in pairs(beta_values)
    println(" Working on T = ", temp[j])
    P_MH_uniform = MH_uniform(N,H, beta)
    e,v  = eigs(P_MH_uniform, nev = 2, which=:LM)
    gap_MH[j] = abs(1-abs(e[2]))
    P_MH_loc = MH_local(N,H, beta)
    e,v  = eigs(P_MH_loc, nev = 2, which=:LM)
    gap_MH_loc[j] = abs(1-abs(e[2]))
end
save_object("Data/Ising-Chain/Classical/N"*string(N)*"MHUniformTemp", gap_MH)
save_object("Data/Ising-Chain/Classical/N"*string(N)*"MHLocTemp", gap_MH_loc)


# get temperature data
# for (i,beta) in pairs(beta_values)
#     println("Working on beta ", beta)
#     h = 0
#     couplings = ones(N)
#     # couplings[end] = 0
#     H = ising_ham(N, couplings, h)
#     num_values = 10
#     kappa_values = range(0,50, length=num_values)
#     eta_values = range(0,50, length=num_values)
#     gap_all = zeros(num_values,num_values)
#     for kappa_i in (1:num_values)
#         println("kappa: ", kappa_i)
#         for eta_i in (1:num_values)
#             println("  eta: ", eta_i)
#             Hm = mixing_ham(N) 
#             P = mixing_qHMC(N, H,kappa_values[kappa_i],eta_values[eta_i],beta, Hm) 
#             e,v  = eigs(P, nev = 2, which=:LM)
#             gap_all[kappa_i,eta_i]  = abs(1-abs(e[2]))   
#         end
#     end
#     name = "Data/Ising-Chain/qMCMC/Grid-Search/Temperature/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
#     save_object(name, gap_all)
# end

# # Average temperature gap
gap_av = zeros(num_points,2)
gap_best = zeros(num_points)
num_values = 50
for (i,beta) in pairs(beta_values)
    println(beta)
    name = "Data/Ising-Chain/qMCMC/Grid-Search/Temperature/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
    gap_all= load_object(name)
    max,cord = findmax(gap_all)
    gap_best[i] = max
    gap_av[i,1] = mean(gap_all)
    gap_av[i,2] = stdm(gap_all, mean(gap_all))
end

save_object("Data/Ising-Chain/qMCMC/N"*string(N)*"AverageTemp", gap_av)
save_object("Data/Ising-Chain/qMCMC/N"*string(N)*"BestTemp", gap_best)