using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2

N = 9
beta = 5
num_values = 300

name = "Data/Ising-Chain/qMCMC/Grid-Search/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
gap_all= load_object(name)

max_param = 50
# max,cord = findmax(gap_all)
# println((cord[2]-1)/num_values*max_param)
# println((cord[1]-1)/num_values*max_param)


# plt.title(L"qHMC Gap Ising $\beta=$ "*string(beta)*" N = "*string(N))
plt.imshow(gap_all, origin="lower",extent = [0, max_param  , 0, max_param ],aspect="auto")
# plt.clim(0, 0.1) 
bar = plt.colorbar()
# plt.scatter((cord[2]-1)/num_values*max_param, (cord[1]-1)/num_values*max_param, color="red", marker=".")

plt.xlabel(L"$\eta$")
plt.ylabel(L"$\kappa$")
bar.set_label(L"$\delta$")
name = "Figures/Ising-Chain/qMCMC/GridSearch/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".svg"
plt.savefig(name)
plt.show()