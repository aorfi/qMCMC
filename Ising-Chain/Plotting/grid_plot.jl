using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2
using Statistics

N = 8
beta = 1000.0
num_values = 50

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
# name = "Figures/Ising-Chain/qMCMC/GridSearch/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".svg"
# plt.savefig(name)
plt.show()



# name = "Data/Ising-Chain/qMCMC/Grid-Search/Temperature/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
# gap_all= load_object(name)
# log_data = -log.(1 .- gap_all)
# m = mean(log_data)
# st_dev = stdm(log_data,mean(log_data))
# st_err = st_dev/sqrt(length(log_data))
# lb = exp(log(m)-st_dev)
# ub = exp(log(m)+st_dev)

# plt.hist(vcat(gap_all...), bins = 100)
# plt.axvline(m, color= "black", label = "Average ")
# # plt.axvline(m+st_dev, color= "black", label = "Average ")
# # plt.axvline(m-st_dev, color= "black", label = "Average ")
# plt.axvline(lb, color= "red")
# plt.axvline(ub, color= "red")
# plt.xlabel(L"$\delta$")
# plt.legend()
# plt.show()