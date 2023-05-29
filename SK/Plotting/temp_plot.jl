using LaTeXStrings
using PythonCall
plt = pyimport("matplotlib.pyplot")
using JLD2
using LsqFit
cd(pwd, "../..")


N = 8
# # plot gap vs temp for fixed N 
temp = 10 .^ (range(-3,stop=3,length=50))
gapMH = load_object("Data/SK/Classical/N"*string(N)*"MHTempAv")[:,1]
gapMH_std = load_object("Data/SK/Classical/N"*string(N)*"MHTempAv")[:,2]
gapMHl = load_object("Data/SK/Classical/N"*string(N)*"MHLocTempAv")[:,1]
gapMHl_std = load_object("Data/SK/Classical/N"*string(N)*"MHLocTempAv")[:,2]
gap_av = load_object("Data/SK/qMCMC/N"*string(N)*"AverageTemp")[:,1]
gap_std = load_object("Data/SK/qMCMC/N"*string(N)*"AverageTemp")[:,2]
gap_best = load_object("Data/SK/qMCMC/N"*string(N)*"BestTemp")[:,1]
gap_best_std = load_object("Data/SK/qMCMC/N"*string(N)*"BestTemp")[:,2]


plt.plot(temp, gapMH, label = "Uniform Proposal", color = "tab:green")
plt.fill_between(temp, gapMH-gapMH_std, gapMH+gapMH_std, color = "tab:green", alpha = 0.2)
plt.plot(temp, gapMHl, label = "Local Proposal", color = "tab:blue")
plt.fill_between(temp, gapMHl-gapMHl_std, gapMHl+gapMHl_std, color = "tab:blue", alpha = 0.2)
plt.plot(temp, gap_av, label = "Average Quantum Proposal", color = "tab:purple")
plt.fill_between(temp, gap_av-gap_std , gap_av+gap_std , color = "tab:purple", alpha = 0.2)
plt.plot(temp, gap_best, label = "Best Quantum Proposal", color = "tab:red")
plt.fill_between(temp, gap_best-gap_best_std , gap_best+gap_best_std , color = "tab:red", alpha = 0.2)

# plt.title( "SK Fully Connected N = "*string(N))
plt.ylabel(L"Average $\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.ylim(10^(-4),2)
plt.grid("both","major")
plt.legend()
name = "Figures/SK/TempN"*string(N)*".svg"
plt.savefig(name)
plt.show()



