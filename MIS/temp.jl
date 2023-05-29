using Arpack
using JLD2
using Random, Distributions
using Statistics

include("load_data.jl")
include("MIS_classical_P.jl")
include("../classical_P.jl")
include("../quantum_P.jl")
include("MIS_hamiltonian.jl")

a = 4.5
Rb = 7.5

side_len = 4
D = 0.5
num_points = 50
temp = 10 .^ (range(-3,stop=6,length=num_points))
beta_values = 1 ./ temp
configurations = load_object("Data/MIS/Configurations/L"*string(side_len)*"D"*string(D)*"config")
runs = 100


# gap_MH_all = zeros(num_points,2)
# gap_MH_loc_all = zeros(num_points,2)
# gap_RSA_all = zeros(num_points,2)
# gap_SA_all = zeros(num_points,2)

# for (i,beta) in pairs(beta_values)
#     println(" Working on T = ", temp[i])
#     gap_MH = zeros(runs)
#     gap_MH_loc = zeros(runs)
#     gap_RSA = zeros(runs)
#     gap_SA = zeros(runs)

#     for j in (1:runs)
#         println(" Working on run = ",j)
#         graph_mask = configurations[j]["graph_mask"]
#         atoms = generate_sites(SquareLattice(), side_len, side_len; scale =a)
#         atoms = apply_graph_mask(atoms, graph_mask) 
#         N = length(atoms)
#         ν = 100
#         H = MIS_ham(atoms, ν)
        
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

#         P_ryd = P_rydberg_SA(H, atoms, beta, Rb)
#         try
#             e,v  = eigs(P_ryd , nev = 2, which=:LM)
#             gap_RSA[j] = abs(1-abs(e[2]))
#         catch
#             e,v = eigen(Matrix(P_ryd))
#             gap_RSA[j]= abs(1-abs(e[end-1]))
#         end

#         ϵ = 0.2
#         P_MIS = P_MIS_SA(H, ϵ,atoms, beta, Rb)
#         try
#             e,v  = eigs(P_MIS, nev = 2, which=:LM)
#             gap_SA[j] = abs(1-abs(e[2]))
#         catch
#             e,v = eigen(Matrix(P_MIS))
#             gap_SA[j]= abs(1-abs(e[end-1]))
#         end
#     end
#     gap_MH_all[i,1] = mean(gap_MH)
#     gap_MH_all[i,2] = stdm(gap_MH, mean(gap_MH))
#     gap_MH_loc_all[i,1] = mean(gap_MH_loc)
#     gap_MH_loc_all[i,2] = stdm(gap_MH_loc, mean(gap_MH_loc))
#     gap_RSA_all[i,1] = mean(gap_RSA)
#     gap_RSA_all[i,2] = stdm(gap_RSA, mean(gap_RSA))
#     gap_SA_all[i,1] = mean(gap_SA)
#     gap_SA_all[i,2] = stdm(gap_SA, mean(gap_SA))
# end
# save_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"MHTempAv", gap_MH_all)
# save_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"MHLocTempAv", gap_MH_loc_all)
# save_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"RSATempAv", gap_RSA_all)
# save_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"SATempAv", gap_SA_all)


# get quantum data

for (i,beta) in pairs(beta_values)
    println("Working on beta ", beta)
    num_values = 10
    kappa_values = range(0,50, length=num_values)
    eta_values = range(0,50, length=num_values)
    gap_all = zeros(num_values,num_values,runs)
    gap_all_ryd= zeros(num_values,num_values,runs)
    for kappa_i in (1:num_values)
        println("kappa: ", kappa_i)
        for eta_i in (1:num_values)
            println("  eta: ", eta_i)
            Hm = mixing_ham(N) 
            for j in (1:runs)
                # println(" Working on run = ",j)
                graph_mask = configurations[j]["graph_mask"]
                atoms = generate_sites(SquareLattice(), side_len, side_len; scale =a)
                atoms = apply_graph_mask(atoms, graph_mask) 
                N = length(atoms)
                ν = 100
                H = MIS_ham(atoms, ν)

                Hc = MIS_ham(atoms, ν)
                Hm = mixing_ham(N)
                P_qHMC = mixing_qHMC(N, Hc,kappa_values[kappa_i],eta_values[eta_i],beta,Hm)
                try
                    e,v  = eigs(P_qHMC, nev = 2, which=:LM)
                    gap_all[kappa_i,eta_i,j] = 1-abs(e[2])
                catch
                    println("Arpack method out of iteration")
                    e,v  = eigen(Matrix(P_qHMC))
                    gap_all[kappa_i,eta_i,j] = 1-abs(e[end-1])
                end

                Hm = mixing_exchange(atoms, Rb)
                P_qHMC_ryd = mixing_rydberg_SA_qHMC(Hc, atoms,kappa_values[kappa_i],eta_values[eta_i],beta,Hm, Rb)
                try
                    e,v  = eigs(P_qHMC_ryd, nev = 2, which=:LM)
                    gap_all_ryd[kappa_i,eta_i,j] = 1-abs(e[2])
                catch
                    println("Arpack method out of iteration")
                    e,v  = eigen(Matrix(P_qHMC_ryd))
                    gap_all_ryd[kappa_i,eta_i,j] = 1-abs(e[end-1])
                end
            end
        end
    end
    save_object("Data/MIS/qMCMC/Grid-Search/Temperature/"*string(num_values)*"L"*string(side_len)*"D"*string(D)*"beta"*string(beta), gap_all)
    save_object("Data/MIS/qMCMC/Grid-Search/Temperature/RSA/"*string(num_values)*"L"*string(side_len)*"D"*string(D)*"beta"*string(beta), gap_all_ryd)
end

