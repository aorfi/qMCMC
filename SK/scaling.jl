using Arpack
using JLD2
using Random, Distributions
using Statistics
Random.seed!(123)

include("SK_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")


N_values = (4:10) 
num = length(N_values)
beta = 5
runs = 100


gap_all_MH = zeros(num,2)
gap_all_MH_loc = zeros(num,2)

for i in (1:num)
    N = N_values[i]
    println(" Working on N = ",N)
    gap_MH = zeros(runs)
    gap_MH_loc = zeros(runs)
    for j in (1:runs)
        println(" Working on run = ",j)
        h = rand(Normal(0,1),N)
        couplings = rand(Normal(0,1),sum(1:N-1))
        H = SK_ham(N,couplings,h)
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
    end
    gap_all_MH[i,1] = mean(gap_MH)
    gap_all_MH[i,2] = stdm(gap_MH, mean(gap_MH))
    gap_all_MH_loc[i,1] = mean(gap_MH_loc)
    gap_all_MH_loc[i,2] = stdm(gap_MH_loc, mean(gap_MH_loc))
end
save_object("Data/SK/Classical/MHScaling", gap_all_MH)
save_object("Data/SK/Classical/MHLocScaling", gap_all_MH_loc)




