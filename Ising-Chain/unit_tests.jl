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

N = 4
h = 0
couplings = ones(N)
# couplings[end] = 0
H = ising_ham(N, couplings, h)
# display(Matrix(H))
energies = [H[i,i] for i in (1:2^N)]
println("energies: ", energies)
# display(counter(energies))
beta = 1
Z = sum([exp(-beta*e) for e in energies])
println("Z ", Z)
# println(2*(2*cosh(beta))^(N-1))
# println((2cosh(beta))^N+(2sinh(beta))^N)


# # beta = 1
# # gs = normalize(exp.(-(beta*energies)))
# gs = zeros(2^N)
# gs[1] = 1
# # gs[end] = 1
# indicies = []
# mixed = mixing_ham(N)*gs
# for i in (1:2^N)
#     if mixed[i] ==1
#         append!(indicies,i)
#     end
# end
# mixed_energies = H*(mixing_ham(N)*gs)
# # println(indicies)
# # println(mixing_ham(N)*gs)
# # println(H*(mixing_ham(N)*gs))
# for i in indicies
#     println(mixed_energies[i])
# end







# # Check construction of P matrices

# P_glab = glab_uniform(N,H, beta)
# eg,vg = eigen(Matrix(transpose(P_glab )))
# # overlap of steady state
# println("Overlap Glauber: ",dot(gs,vg[:,end]))
# gapg = abs(1-abs(eg[end-1]))
# println("Gap glaub: ", gapg)

# P_MH = MH_local(N,H, beta)
# eMH,vMH = eigen(Matrix(transpose(P_MH)))
# println("overlap MH: ",dot(gs,vMH[:,end]))
# gapMH = abs(1-abs(eMH[end-1]))
# println("Gap MH: ", gapMH)

# α = 2
# η = 3
# Hm = mixing_ham(N) 
# P_qHMC = mixing_qHMC(N, H,α,η,beta, Hm)
# eqHMC,vqHMC  = eigen(Matrix(transpose(P_qHMC)))
# eq,vq = eigen(Matrix(P_qHMC))
# println("overlap qHMC: ",dot(gs,vqHMC[:,end]))
# gap = abs(1-abs(eqHMC[end-1]))
# println("Gap qHMC: ", gap)

