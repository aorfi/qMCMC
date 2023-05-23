using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using Random, Distributions

include("SK_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")
include("../conductance.jl")

beta = 5
N = 3
h = rand(Normal(0,1),N)
couplings = rand(Normal(0,1),sum(1:N-1))
H = SK_ham(N,couplings,h)
energies = [H[i,i] for i in (1:2^N)]
println("energies: ", energies)
Z = sum([exp(-beta*e) for e in energies])
gs = exp.(-(beta*energies))/Z

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

κ = 2
η = 3
Hm = mixing_ham(N) 
P_qHMC = mixing_qHMC(N, H,κ,η,beta, Hm)
println("gs*P - gs: ",transpose(gs)*P_qHMC - transpose(gs))
eqHMC,vqHMC  = eigen(Matrix(transpose(P_qHMC)))
gap = abs(1-abs(eqHMC[end-1]))
println("Gap qHMC: ", gap)
