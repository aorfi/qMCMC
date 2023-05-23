using Arpack
using JLD2
using Random, Distributions
using Statistics
Random.seed!(123)

include("SK_hamiltonian.jl")
include("../classical_P.jl")
include("../quantum_P.jl")


N = 8
runs = 100
num_points = 50
temp = 10 .^ (range(-3,stop=3,length=num_points))
beta_values = 1 ./ temp


gap_all_MH = zeros(num_points,2)
gap_all_MH_loc = zeros(num_points,2)

for (i,beta) in pairs(beta_values)
    println(" Working on T = ", temp[i])
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
save_object("Data/SK/Classical/N"*string(N)*"MHTempAv", gap_all_MH)
save_object("Data/SK/Classical/N"*string(N)*"MHLocTempAv", gap_all_MH_loc)



