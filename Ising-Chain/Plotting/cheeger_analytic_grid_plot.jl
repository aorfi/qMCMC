using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2

include("../Cheeger/cheeger_analytic.jl")


N = 6
uniform_gap = 2.0^(-N+1)
println("Uniform Gap ", uniform_gap )
beta = 5
num_values = 300
kappa_values = range(0.00000001,50, length=num_values)
eta_values = range(0.00000001,50, length=num_values)
phi = zeros(num_values,num_values)
phi_lh =[]
for kappa_i in (1:num_values)
    # println("kappa: ", kappa_i)
    for eta_i in (1:num_values)
        # println("  eta: ", eta_i)
        phi[kappa_i,eta_i]= all_but_ground_states(N,kappa_values[kappa_i],eta_values[eta_i])
        # phi_2 = one_ground_state(N,kappa_values[kappa_i],eta_values[eta_i])
        # phi[kappa_i,eta_i] = min(phi_1,phi_2)
        # if eta_values[eta_i]>kappa_values[kappa_i]
        #     append!(phi_lh, min(phi_1,phi_2) )
        # end
    end
end
println( "average phi ",mean(phi))
# println( "average phi lower half ",mean(phi_lh))



max_param = 50
plt.imshow(phi, origin="lower",extent = [0, max_param  , 0, max_param ],aspect="auto")
# plt.clim(0, uniform_gap) 
bar = plt.colorbar()
plt.xlabel(L"$\eta$")
plt.ylabel(L"$\kappa$")
bar.set_label(L"$\delta$")
# name = "Figures/Ising-Chain/qMCMC/GridSearch/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".svg"
# plt.savefig(name)
plt.show()