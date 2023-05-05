using JLD2

include("ising_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")

N = 14
h = 0
couplings = ones(N)
# couplings[end] = 0
H = ising_ham(N, couplings, h)

num_points = 50
temp = 10 .^ (range(-3,stop=3,length=num_points))
beta_values = 1 ./ temp

# gap_glab = zeros(num_points)
# gap_glab_loc = zeros(num_points)
gap_MH = zeros(num_points)
gap_MH_loc = zeros(num_points)
# gap_qHMC = zeros(num_points)
# gap_rand = zeros(num_points)
for (j,beta) in pairs(beta_values)
    println(" Working on T = ", temp[j])
    # M_glab = mixing_glab(N,H, beta)
    # try
    #     e,v  = eigs(M_glab, nev = 2, which=:LM)
    #     gap_glab[j] = abs(1-abs(e[2]))
    # catch
    #     e,v = eigen(Matrix(M_glab))
    #     gap_glab[j]= abs(1-abs(e[end-1]))
    # end
    # M_glab_loc = mixing_glab_loc(N,H, beta)
    # try
    #     e,v  = eigs(M_glab_loc, nev = 2, which=:LM)
    #     gap_glab_loc[j] = abs(1-abs(e[2]))
    # catch
    #     e,v = eigen(Matrix(M_glab_loc))
    #     gap_glab_loc[j]= abs(1-abs(e[end-1]))
    # end
    P_MH_uniform = MH_uniform(N,H, beta)
    e,v  = eigs(P_MH_uniform, nev = 2, which=:LM)
    gap_MH[j] = abs(1-abs(e[2]))

    P_MH_loc = MH_local(N,H, beta)
    e,v  = eigs(P_MH_loc, nev = 2, which=:LM)
    gap_MH_loc[j] = abs(1-abs(e[2]))

    # Hm = mixing_ham(N) 
    # M_qHMC = mixing_qHMC(N, H,α,η,beta, Hm)
    # try
    #     e,v  = eigs(M_qHMC, nev = 2, which=:LM)
    #     gap_qHMC[j] = abs(1-abs(e[2]))
    # catch
    #     e,v = eigen(Matrix(M_qHMC ))
    #     gap_qHMC[j]= abs(1-abs(e[end-1]))
    # end
    # M_rand = random_mixing(N,H,beta,Hm, 100)
    # try
    #     e,v  = eigs(M_rand, nev = 2, which=:LM)
    #     gap_rand[j] = abs(1-abs(e[2]))
    # catch
    #     e,v = eigen(Matrix(M_rand))
    #     gap_rand[j]= abs(1-abs(e[end-1]))
    # end
end
# save_object("Data/Ising/Classical/N"*string(N)*"GlaubTemp", gap_glab)
# save_object("Data/Ising/Classical/N"*string(N)*"GlaubLocTemp", gap_glab_loc)
save_object("Data/Ising-Chain/Classical/N"*string(N)*"MHUniformTemp", gap_MH)
save_object("Data/Ising-Chain/Classical/N"*string(N)*"MHLocTemp", gap_MH_loc)
# save_object("Data/Ising/qHMC/N"*string(N)*"TempBest", gap_qHMC)
# save_object("Data/Ising/qHMC/RandN"*string(N)*"Temp", gap_rand)

