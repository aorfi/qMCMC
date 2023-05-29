using LaTeXStrings
using PythonCall
plt = pyimport("matplotlib.pyplot")
using JLD2
using LsqFit
cd(pwd, "../..")


side_len = 4
D = 0.5
# # plot gap vs temp for fixed N 
temp = 10 .^ (range(-3,stop=6,length=50))
gapMH = load_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"MHTempAv")[:,1]
gapMH_std = load_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"MHTempAv")[:,2]
gapMHl = load_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"MHLocTempAv")[:,1]
gapMHl_std = load_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"MHLocTempAv")[:,2]
gapRSA = load_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"RSATempAv")[:,1]
gapRSA_std = load_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"RSATempAv")[:,2]
gapSA = load_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"SATempAv")[:,1]
gapSA_std = load_object("Data/MIS/Classical/L"*string(side_len)*"D"*string(D)*"SATempAv")[:,2]
println(gapRSA)

plt.plot(temp, gapMH, label = "Uniform Proposal", color = "tab:green")
plt.fill_between(temp, gapMH-gapMH_std, gapMH+gapMH_std, color = "tab:green", alpha = 0.2)
plt.plot(temp, gapMHl, label = "Local Proposal", color = "tab:blue")
plt.fill_between(temp, gapMHl-gapMHl_std, gapMHl+gapMHl_std, color = "tab:blue", alpha = 0.2)
plt.plot(temp, gapRSA, label = "Exchange Proposal", color = "black")
plt.fill_between(temp, gapRSA-gapRSA_std, gapRSA+gapRSA_std, color = "black", alpha = 0.2)


plt.ylabel(L"Average $\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.ylim(10^(-3),2)
plt.grid("both","major")
plt.legend()
name = "Figures/MIS/TempL"*string(side_len)*"D"*string(D)*".svg"
plt.savefig(name)
plt.show()



