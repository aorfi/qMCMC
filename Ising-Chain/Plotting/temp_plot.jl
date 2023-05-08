using LaTeXStrings
using PythonCall
plt = pyimport("matplotlib.pyplot")
using JLD2
cd(pwd, "../..")


N = 9
# # plot gap vs temp for fixed N 
temp = 10 .^ (range(-3,stop=3,length=50))
gapMH = load_object("Data/Ising-Chain/Classical/N"*string(N)*"MHUniformTemp")
gapMHl = load_object("Data/Ising-Chain/Classical/N"*string(N)*"MHLocTemp")
# gapG = load_object("Data/Ising/Classical/N"*string(N)*"GlaubTemp")
# gapGl = load_object("Data/Ising/Classical/N"*string(N)*"GlaubLocTemp")
# gapq = load_object("Data/Ising/qHMC/N"*string(N)*"TempBest")
# gapRand = load_object("Data/Ising/qHMC/RandN"*string(N)*"Temp")


# plt.scatter(temp, gapG, label = "Glaubler Uniform", color = "tab:blue")
# plt.scatter(temp, gapGl, label = "Glaubler Local", color = "tab:orange")
plt.scatter(temp, gapMH, label = "MH Uniform", color = "tab:green")
plt.scatter(temp, gapMHl, label = "MH Local", color = "tab:red")
# plt.scatter(temp, gapq, label = "qHMC Best", color = "tab:purple")
# plt.scatter(temp, gapRand, label = "qHMC Rand", color = "tab:brown")


plt.title( "1D Ising OBC N = "*string(N))
plt.ylabel(L" $\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.ylim(10^(-5),1.5)
plt.grid("both","major")
plt.axvline(1/5)
plt.legend()
# plt.savefig("Figures/Ising/qHMC/GapTempN"*string(N)*"beta"*string(beta)*".png")
plt.show()



