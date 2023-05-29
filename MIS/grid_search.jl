using LinearAlgebra
using Bloqade
using PythonCall
using Arpack
using Graphs
using JLD2
using GenericTensorNetworks


include("load_data.jl")
include("MIS_classical_P.jl")
include("MIS_quantum_P.jl")
include("../classical_P.jl")
include("../quantum_P.jl")
include("MIS_hamiltonian.jl")
cd(pwd, "..")

a = 4.5
Rb = 7.5

runs = 1
num = 6
N_all = [4,5,6,7,8,10,11,13,15,18,20]
L = [3,3,3,3,4,4,4,4,5,5,5]
D_all = [0.6,0.4,0.3,0.2,0.5,0.4,0.3,0.2,0.4,0.3,0.2]
config_all = Matrix{Any}(undef, num, runs)
for i in (1:num)
    for j in (1:runs)
        config_all[i,j] = load_object("Data/MIS/Configurations/L"*string(L[i])*"D"*string(D_all[i])*"config")[j]
    end 
end

for i in (1:1)
    beta = 10
    num_values = 20
    eta_values = range(0,50, length=num_values)
    kappa_values = range(0,50, length=num_values)
    gap_all = zeros(num_values,num_values,runs)
    gap_all_ryd = zeros(num_values,num_values,runs)
    for kappa_i in (1:num_values)
        for eta_i in (1:num_values)
            println("  kappa: ", kappa_i)
            println("  eta: ", eta_i)
            for j in (1:runs)
                # println(" Working on run = ",j)
                graph_mask = config_all[i,j]["graph_mask"]
                atoms = generate_sites(SquareLattice(), L[i], L[i]; scale =a)
                atoms = apply_graph_mask(atoms, graph_mask) 
                N = length(atoms)
    
                ν = 100
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
    name = "Data/MIS/qMCMC/Grid-Search/"*string(num_values)*"L"*string( L[i])*"D"*string(D_all[i])*"beta"*string(beta)*"Av"*string(runs)
    name_ryd = "Data/MIS/qMCMC/Grid-Search/RSA/"*string(num_values)*"L"*string( L[i])*"D"*string(D_all[i])*"beta"*string(beta)*"Av"*string(runs)
    save_object(name, gap_all)
    save_object(name_ryd, gap_all_ryd)
end



