using Statistics

include("../ising_hamiltonian.jl")
include("../../classical_P.jl")
include("../../conductance.jl")
include("../../quantum_P.jl")


function one_ground_state(N,κ,η)
    k_pos0 = (2*pi/N)*range(1,(N)/2) .- pi/N
    U0_anal_list = []
    for k in k_pos0
        ek = sqrt((η-κ*cos(k))^2+(κ*sin(k))^2)
        theta_k = acos((η-κ*cos(k))/ek)
        append!(U0_anal_list,cos(0.5*(k+theta_k))^2*exp(-2im*ek)+sin(0.5*(k+theta_k))^2*exp(2im*ek))
    end
    U0_anal = prod(U0_anal_list)
    k_pos1 = (2*pi/N)*range(1,(N)/2-1) 
    U1_anal_list = []
    for k in k_pos1
        ek = sqrt((η-κ*cos(k))^2+(κ*sin(k))^2)
        theta_k = acos((η-κ*cos(k))/ek)
        append!(U1_anal_list,cos(0.5*(k+theta_k))^2*exp(-2im*ek)+sin(0.5*(k+theta_k))^2*exp(2im*ek))
    end
    U1_anal = prod(U1_anal_list*exp(1im*2*κ))
    U00 = U0_anal*conj(U0_anal)
    U11 = U1_anal*conj(U1_anal)
    U01 = U0_anal*conj(U1_anal)
    U10 = U1_anal*conj(U0_anal)
    phi_anal  = 0.25*(U00+U11-U01-U10)
    return real(phi_anal)
end

function all_but_ground_states(N,κ,η)
    k_pos0 = (2*pi/N)*range(1,(N)/2) .- pi/N
    overlapU_k = []
    overlapU_others = []
    for k in k_pos0
        ek = sqrt((η-κ*cos(k))^2+(κ*sin(k))^2)
        theta_k = acos((η-κ*cos(k))/ek)
        append!(overlapU_others , 1- η^2*sin(k)^2*sin(2*ek)^2/ek^2 )
        append!(overlapU_k , η^2*sin(2*ek)^2*sin(k)^2/ek^2 )
    end
    overlapU_anal = []
    for i in (1:length(overlapU_k))
        append!(overlapU_anal, overlapU_k[i]*prod(overlapU_others)/overlapU_others[i] )
    end
    p0_over = sum(overlapU_anal)

    k_pos1 = (2*pi/N)*range(1,(N)/2-1)
    overlapU_k = []
    overlapU_others = []
    for k in k_pos1
        ek = sqrt((η-κ*cos(k))^2+(κ*sin(k))^2)
        theta_k = acos((η-κ*cos(k))/ek)
        append!(overlapU_others , 1- η^2*sin(k)^2*sin(2*ek)^2/ek^2 )
        append!(overlapU_k , η^2*sin(2*ek)^2*sin(k)^2/ek^2 )
    end
    overlapU_anal = []
    for i in (1:length(overlapU_k))
        append!(overlapU_anal, overlapU_k[i]*prod(overlapU_others)/overlapU_others[i] )
    end
    p1_over = sum(overlapU_anal)
    phi_anal = (p0_over+p1_over)/(N*(N-1))
    return real(phi_anal)
end


#Unit tests
# beta = 5
# N = 4
# κ = 4
# η = 3
# couplings = ones(N)
# Hc = ising_ham(N, couplings, 0)
# Hm = mixing_ham(N) 
# P = mixing_qHMC(N, Hc,κ,η,beta, Hm) 
# A = (2:2^N-1)
# phi_A = conductance(Hc, P, N, A,beta)
# println("Phi ", phi_A)
# println("Phi analytic ", all_but_ground_states(N,κ,η))
# A = [1]
# phi_A = conductance(Hc, P, N, A,beta)
# println("Phi ", phi_A)
# println("Phi analytic ", one_ground_state(N,κ,η))


# N_values = (50:4:150)
# phi_N = zeros(length(N_values),2)
# for i in (1:length(N_values))
#     N = N_values[i]
#     println("N: ", N)
#     beta = 5
#     num_values = 100
#     # kappa_values = range(10,50, length=num_values)
#     # eta_values = range(10,50, length=num_values)
#     kappa_values = [10]
#     eta_values =[40]
#     phi = zeros(length(kappa_values),length(eta_values))
#     for kappa_i in (1:length(kappa_values))
#         # println("kappa: ", kappa_i)
#         for eta_i in (1:length(eta_values))
#             # println("  eta: ", eta_i)
#             # if eta_values[eta_i]>(kappa_values[kappa_i])
#                 phi_1 = all_but_ground_states(N,kappa_values[kappa_i],eta_values[eta_i])
#                 phi_2 = one_ground_state(N,kappa_values[kappa_i],eta_values[eta_i])
#                 phi[kappa_i,eta_i] = min(phi_1,phi_2)
#             # end
#         end
#     end
#     phi_N[i,1] = mean(phi)
#     phi_N[i,2] = stdm(phi, mean(phi))
# end
# println(phi_N[:,1])
# # name = "Data/Ising-Chain/qMCMC/Cheeger/AverageBottleneckScaling"
# # name = "Data/Ising-Chain/qMCMC/Cheeger/LowerHalfBottleneckScaling"
# name = "Data/Ising-Chain/qMCMC/Cheeger/BottleneckScalingEta40Kappa20"
# save_object(name, phi_N)
