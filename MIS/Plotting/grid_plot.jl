using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2
using Statistics

N = 7
beta = 10
num_values = 20
runs = 1
index = 1
N_all = [4,5,6,7,8,10,11,13,15,18,20]
L = [3,3,3,3,4,4,4,4,5,5,5]
D_all = [0.6,0.4,0.3,0.2,0.5,0.4,0.3,0.2,0.4,0.3,0.2]

# name = "Data/MIS/qMCMC/Grid-Search/"*string(num_values)*"L"*string( L[index])*"D"*string(D_all[index])*"beta"*string(beta)*"Av"*string(runs)
name = "Data/MIS/qMCMC/Grid-Search/RSA/"*string(num_values)*"L"*string( L[index])*"D"*string(D_all[index])*"beta"*string(beta)*"Av"*string(runs)

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
# plt.imshow(gap_all[:,:,3], origin="lower",extent = [0, max_param  , 0, max_param ],aspect="auto")
# plt.clim(0, 0.1) 
bar = plt.colorbar()
# plt.scatter((cord[2]-1)/num_values*max_param, (cord[1]-1)/num_values*max_param, color="red", marker=".")

plt.xlabel(L"$\eta$")
plt.ylabel(L"$\kappa$")
bar.set_label(L"Average $\delta$")
# name = "Figures/SK/qMCMC/GridSearch/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".svg"
# plt.savefig(name)
plt.show()