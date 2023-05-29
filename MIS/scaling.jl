using LinearAlgebra
using Bloqade
using PythonCall
using Arpack
using Graphs
using JLD2
using GenericTensorNetworks


include("load_data.jl")
include("MIS_classical_P.jl")
include("../classical_P.jl")
include("../quantum_P.jl")
include("MIS_hamiltonian.jl")
cd(pwd, "..")

a = 4.5
Rb = 7.5
beta = 15

# # Scaling averaged over configurations
# runs = 100
# num = 6
# N_all = [4,5,6,7,8,10,11,13,15,18,20]
# L = [3,3,3,3,4,4,4,4,5,5,5]
# D_all = [0.6,0.4,0.3,0.2,0.5,0.4,0.3,0.2,0.4,0.3,0.2]
# config_all = Matrix{Any}(undef, num, runs)
# for i in (1:num)
#     for j in (1:runs)
#         config_all[i,j] = load_object("Data/MIS/Configurations/L"*string(L[i])*"D"*string(D_all[i])*"config")[j]
#     end 
# end

# gap_MH_all = zeros(num_points,2)
# gap_MH_loc_all = zeros(num_points,2)
# gap_RSA_all = zeros(num_points,2)
# gap_SA_all = zeros(num_points,2)
# for i in (1:num)
#     println("working on N ", N_all[i])
#     gap_MH = zeros(runs)
#     gap_MH_loc = zeros(runs)
#     gap_RSA = zeros(runs)
#     gap_SA = zeros(runs)
#     for j in (1:runs)
#         println(" Working on run = ",j)
#         graph_mask = config_all[i,j]["graph_mask"]
#         atoms = generate_sites(SquareLattice(), L[i], L[i]; scale =a)
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
# save_object("Data/MIS/Classical/scalingMHBeta"*string(beta)*"Av", gap_MH_all)
# save_object("Data/MIS/Classical/scalingMHLocBeta"*string(beta)*"Av", gap_MH_loc_all)
# save_object("Data/MIS/Classical/scalingRydbergRSABeta"*string(beta)*"Av", gap_RSA_all)
# save_object("Data/MIS/Classical/scalingRydbergSABeta"*string(beta)*"Av", gap_SA_all)



# Scaling for fixed HP and deg
configs = load_object("Data/MIS/Configurations/Deg1HP2configs")
index = 4
N_all = [4,5,6,8,15,18,20][1:index]
L = [3,3,3,4,5,5,5][1:index]

gap_MH = zeros(index)
gap_MH_loc = zeros(index)
gap_RSA = zeros(index)
gap_SA = zeros(index)
for j in (1:index)

    graph_mask = configs[j]["graph_mask"]
    atoms = generate_sites(SquareLattice(), L[j], L[j]; scale =a)
    atoms = apply_graph_mask(atoms, graph_mask) 
    N = length(atoms)
    println("working on N ", N)
    ν = 100
    H = MIS_ham(atoms, ν)
    P_MH_uniform = MH_uniform(N,H, beta)
    try
        e,v  = eigs(P_MH_uniform, nev = 2, which=:LM)
        gap_MH[j] = abs(1-abs(e[2]))
    catch
        e,v = eigen(Matrix(P_MH_uniform))
        gap_MH[j]= abs(1-abs(e[end-1]))
    end

    P_MH_loc = MH_local(N,H, beta)
    try
        e,v  = eigs(P_MH_loc, nev = 2, which=:LM)
        gap_MH_loc[j] = abs(1-abs(e[2]))
    catch
        e,v = eigen(Matrix(P_MH_loc))
        gap_MH_loc[j]= abs(1-abs(e[end-1]))
    end

    P_ryd = P_rydberg_SA(H, atoms, beta, Rb)
    try
        e,v  = eigs(P_ryd , nev = 2, which=:LM)
        gap_RSA[j] = abs(1-abs(e[2]))
    catch
        e,v = eigen(Matrix(P_ryd))
        gap_RSA[j]= abs(1-abs(e[end-1]))
    end

    ϵ = 0.2
    P_MIS = P_MIS_SA(H, ϵ,atoms, beta, Rb)
    try
        e,v  = eigs(P_MIS, nev = 2, which=:LM)
        gap_SA[j] = abs(1-abs(e[2]))
    catch
        e,v = eigen(Matrix(P_MIS))
        gap_SA[j]= abs(1-abs(e[end-1]))
    end
end
save_object("Data/MIS/Classical/scalingFixedMHBeta"*string(beta), gap_MH)
save_object("Data/MIS/Classical/scalingFixedMHLocBeta"*string(beta), gap_MH_loc)
save_object("Data/MIS/Classical/scalingFixedRydbergRSABeta"*string(beta), gap_RSA)
save_object("Data/MIS/Classical/scalingFixedRydbergSABeta"*string(beta), gap_SA)

