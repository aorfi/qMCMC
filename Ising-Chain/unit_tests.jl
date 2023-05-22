using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using Random, Distributions
using PythonCall
using DataStructures
plt = pyimport("matplotlib.pyplot")

include("ising_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")
include("../conductance.jl")

beta = 5
N = 2
h = 0
couplings = ones(N)
# couplings[end] = 0
H = ising_ham(N, couplings, h)
# display(Matrix(H))
energies = [H[i,i] for i in (1:2^N)]
println("energies: ", energies)
# display(counter(energies))

Z = sum([exp(-beta*e) for e in energies])
gs = exp.(-(beta*energies))/Z
# println(gs)
# println(normalize(exp.(-(beta*energies))))






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


ratio = norm(real.(Matrix(Hm)),2)/norm(real.(Matrix(H)), 2)
H_a = κ*ratio*H+η*Hm
display(Matrix(exp(-1im*Matrix(H_a))))