using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2

N = 4
beta = 5
num_values = 300

name = "Data/Ising-Chain/qMCMC/Cheeger/CheegerGrid"*string(num_values)*"N"*string(N)*"beta"*string(beta)
gap_all= load_object(name)
name = "Data/Ising-Chain/qMCMC/Grid-Search/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
phi_all= load_object(name)
max,cord = findmax(2*phi_all)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
im1 = ax1.imshow(gap_all , origin="lower",extent = [0, max_param  , 0, max_param ],aspect="auto")
im1.set_clim(0, max) 
bar1 = plt.colorbar(im1, ax = ax1)


im2 = ax2.imshow(2*phi_all, origin="lower",extent = [0, max_param  , 0, max_param ],aspect="auto")
# cbar2 = plt.colorbar(im2,ax2)
im2.set_clim(0, max) 
bar2 = plt.colorbar(im2,ax = ax2)
bar1.set_label(L"$\delta$")
bar2.set_label(L"$2\Phi$")
ax1.set_xlabel(L"$\eta$")
ax2.set_xlabel(L"$\eta$")
ax1.set_ylabel(L"$\kappa$")

# plt.tight_layout()

name = "Figures/Ising-Chain/qMCMC/Cheeger/Comparison"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".svg"
plt.savefig(name)
plt.show()