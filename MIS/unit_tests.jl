using LinearAlgebra
using Bloqade
using PythonCall
using Arpack
using Graphs
using JLD2
using GenericTensorNetworks
matplotlib = pyimport("matplotlib")
plt = pyimport("matplotlib.pyplot")
BloqadeLattices.DEFAULT_BACKGROUND_COLOR[] = "#FFFFFF"

include("load_data.jl")
include("MIS_classical_P.jl")
include("MIS_quantum_P.jl")
include("../classical_P.jl")
include("../quantum_P.jl")
include("MIS_hamiltonian.jl")


a = 4.5
Rb = 7.5

N_all = [4,5,6,7,8,10,11,13,15,18,20]
L = [3,3,3,3,4,4,4,4,5,5,5]
D_all = [0.6,0.4,0.3,0.2,0.5,0.4,0.3,0.2,0.4,0.3,0.2]
index = 1
configs = load_object("Data/MIS/Configurations/L"*string(L[index])*"D"*string(D_all[index])*"config")
graph_mask = configs[1]["graph_mask"]

graph_mask = configs[index]["graph_mask"]
atoms = generate_sites(SquareLattice(), L[index], L[index]; scale =a)
atoms = apply_graph_mask(atoms, graph_mask)
# atoms = generate_sites(SquareLattice(), 3, 3; scale = a)|> random_dropout(0.4)
# display(Bloqade.plot(atoms, blockade_radius = Rb))

beta = 30
N = length(atoms)
ν = 100
H = MIS_ham(atoms, ν)
energies = [H[i,i] for i in (1:2^N)]
println("energies: ", energies)
Z = sum([exp(-beta*e) for e in energies])
gs = exp.(-(beta*energies))/Z
println(gs)

# # Get MIS solution vector
# graph = BloqadeMIS.unit_disk_graph(atoms, Rb)
# space = independent_set_subspace(graph)
# problem = IndependentSet(graph)
# bloqade_configs = GenericTensorNetworks.solve(problem, ConfigsMax())[].c
# num_mis = length(bloqade_configs)
# configs = zeros(num_mis)
# [configs[i] = parse(Int,reverse(string(bloqade_configs[i])); base=2) for i in (1:num_mis)]
# MIS_vec = zeros(2^N) #vector in computational basis
# [MIS_vec[i+1] = 1/(num_mis) for i in (0:2^N-1) if i in configs] # normalized to match steady state
# println(MIS_vec)


# # Check construction of P matrices

P_glab = glab_uniform(N,H, beta)
# check steady state
println("gs*P - gs: ",transpose(gs)*P_glab - transpose(gs))
eg,vg = eigen(Matrix(transpose(P_glab)))
gapg = abs(1-abs(eg[end-1]))
println("Gap glaub: ", gapg)

P_MH = MH_local(N,H, beta)
println("gs*P - gs: ",transpose(gs)*P_MH - transpose(gs))
eMH,vMH = eigen(Matrix(transpose(P_MH)))
gapMH = abs(1-abs(eMH[end-1]))
println("Gap MH: ", gapMH)

P_ryd = P_rydberg_SA(H, atoms, beta, Rb)
println("gs*P - gs: ",transpose(gs)*P_ryd - transpose(gs))
eryd,vryd = eigen(Matrix(transpose(P_ryd)))
gap_ryd = abs(1-abs(eryd[end-1]))
println("Gap Ryd: ", gap_ryd)

ϵ = 0.2
P_MIS = P_MIS_SA(H, ϵ,atoms, beta, Rb)
println("gs*P - gs: ",transpose(gs)*P_MIS - transpose(gs))
eMIS,vMIS = eigen(Matrix(transpose(P_MIS)))
gap_MIS = abs(1-abs(eMIS[end-1]))
println("Gap MIS: ", gap_MIS)



κ = 10
η = 20
Hm = mixing_ham(N) 
P_qHMC = mixing_qHMC(N, H,κ,η,beta, Hm)
println("gs*P - gs: ",transpose(gs)*P_qHMC - transpose(gs))
eqHMC,vqHMC  = eigen(Matrix(transpose(P_qHMC)))
gap = abs(1-abs(eqHMC[end-1]))
println("Gap qHMC: ", gap)


Hm = mixing_exchange(atoms, Rb)  
P_qHMC_ryd = mixing_rydberg_SA_qHMC(H, atoms,κ,η, beta, Hm, Rb)
println("gs*P - gs: ",transpose(gs)*P_qHMC_ryd - transpose(gs))
eqHMC_ryd,vqHMC_ryd  = eigen(Matrix(transpose(P_qHMC_ryd)))
gap = abs(1-abs(eqHMC_ryd[end-1]))
println("Gap qHMC_ryd: ", gap)
