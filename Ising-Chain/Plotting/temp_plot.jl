using LaTeXStrings
using PythonCall
plt = pyimport("matplotlib.pyplot")
using JLD2
cd(pwd, "../..")


N = 8
# # plot gap vs temp for fixed N 
temp = 10 .^ (range(-3,stop=3,length=50))
gapMH = load_object("Data/Ising-Chain/Classical/N"*string(N)*"MHUniformTemp")
gapMHl = load_object("Data/Ising-Chain/Classical/N"*string(N)*"MHLocTemp")
gap_av = load_object("Data/Ising-Chain/qMCMC/N"*string(N)*"AverageTemp")[:,1]
std_av = load_object("Data/Ising-Chain/qMCMC/N"*string(N)*"AverageTemp")[:,2]
gap_best = load_object("Data/Ising-Chain/qMCMC/N"*string(N)*"BestTemp")

plt.plot(temp, gapMH, label = "Uniform Proposal", color = "tab:green")
plt.plot(temp, gapMHl, label = "Local Proposal", color = "tab:blue")
plt.plot(temp, gap_av, label = "Random Quantum Proposal", color = "tab:purple")
plt.fill_between(temp, exp.(log.(gap_av) .- std_av), exp.(log.(gap_av) .+ std_av), color = "tab:purple", alpha = 0.2)
plt.plot(temp, gap_best, label = "Best Quantum Proposal", color = "tab:red")

# plt.title( "1D Ising OBC N = "*string(N))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.ylim(10^(-3),1.5)
plt.grid("both","major")
# plt.axvline(1/5)
plt.legend()
name = "Figures/Ising-Chain/TempN"*string(N)*".svg"
plt.savefig(name)
plt.show()



