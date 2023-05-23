using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2
using Statistics

N = 4
beta = 5
num_values = 10

name = "Data/SK/qMCMC/Grid-Search/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
gap_all= load_object(name)
gap_av = zeros(num_values,num_values)
for kappa_i in (1:num_values)
    for eta_i in (1:num_values)
        gap_av[kappa_i,eta_i] = mean(gap_all[kappa_i,eta_i,:])
    end
end

max_param = 50
# max,cord = findmax(gap_all)
# println((cord[2]-1)/num_values*max_param)
# println((cord[1]-1)/num_values*max_param)


# plt.title(L"qHMC Gap Ising $\beta=$ "*string(beta)*" N = "*string(N))
plt.imshow(gap_av, origin="lower",extent = [0, max_param  , 0, max_param ],aspect="auto")
# plt.clim(0, 0.1) 
bar = plt.colorbar()
# plt.scatter((cord[2]-1)/num_values*max_param, (cord[1]-1)/num_values*max_param, color="red", marker=".")

plt.xlabel(L"$\eta$")
plt.ylabel(L"$\kappa$")
bar.set_label(L"$\delta$")
# name = "Figures/Ising-Chain/qMCMC/GridSearch/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".svg"
# plt.savefig(name)
plt.show()