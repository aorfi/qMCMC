using GLM
using DataFrames
using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2
cd(pwd, "..")


# # Average over graphs
N_all = [4,5,6,7,8,10,11,13,15,18,20]
L = [3,3,3,3,4,4,4,4,5,5,5]
D = [0.6,0.4,0.3,0.2,0.5,0.4,0.3,0.2,0.4,0.3,0.2]
beta = 10

index = 6
gap_MH = load_object("Data/MIS/Classical/scalingMHBeta"*string(beta)*"Av")[:,1][1:index]
gapMH_std = load_object("Data/MIS/Classical/scalingMHBeta"*string(beta)*"Av")[:,2][1:index]
gap_MH_loc = load_object("Data/MIS/Classical/scalingMHLocBeta"*string(beta)*"Av")[:,1][1:index]
gap_MH_loc_std = load_object("Data/MIS/Classical/scalingMHLocBeta"*string(beta)*"Av")[:,2][1:index]
gap_RSA = load_object("Data/MIS/Classical/scalingRydbergRSABeta"*string(beta)*"Av")[:,1][1:index]
gap_RSA_std = load_object("Data/MIS/Classical/scalingRydbergRSABeta"*string(beta)*"Av")[:,2][1:index]
gap_SA = load_object("Data/MIS/Classical/scalingRydbergSABeta"*string(beta)*"Av")[:,1][1:index]
gap_SA_std = load_object("Data/MIS/Classical/scalingRydbergSABeta"*string(beta)*"Av")[:,2][1:index]

plt.errorbar(N_all[1:index],  gap_MH, yerr = gapMH_std , fmt="o",color = "tab:green", label = "Uniform Proposal")
plt.errorbar(N_all[1:index],  gap_MH_loc, yerr = gap_MH_loc_std , fmt="o",color = "tab:blue", label = "Local Proposal")
plt.errorbar(N_all[1:index],  gap_RSA, yerr = gap_RSA_std , fmt="o",color = "black", label = "Exchange Proposal")


plt.ylabel(L" Average $\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid(which="major", zorder=-1.0)
plt.legend()
plt.show()







# # Fixed HP and deg
# index = 4
# N_all = [4,5,6,8,15,18,20][1:index]
# beta = 10 
# gap_MH = load_object("Data/MIS/Classical/scalingFixedMHBeta"*string(beta))
# gap_MH_loc = load_object("Data/MIS/Classical/scalingFixedMHLocBeta"*string(beta))
# gap_RSA = load_object("Data/MIS/Classical/scalingFixedRydbergRSABeta"*string(beta))
# gap_SA = load_object("Data/MIS/Classical/scalingFixedRydbergSABeta"*string(beta))
 


# # label_scatterRSA = L"Rydberg SA $2^{-kN}$ k= "*string(scaleRSA)
# plt.scatter(N_all  , gap_RSA, color = "black")
# # plt.plot(x,exp.(model(x,paramRSA)),linestyle = "dashed", label = label_scatterRSA,color = "black")

# # label_scatterSA = L"MIS SA $2^{-kN}$ k= "*string(scaleSA)
# plt.scatter(N_all  , gap_SA, color = "tab:gray")
# # plt.plot(x,exp.(model(x,paramSA)),linestyle = "dashed", label = label_scatterSA,color = "tab:gray")

# # label_scatterMH = L"Uniform MH $2^{-kN}$ k= "*string(scaleMH)
# plt.scatter(N_all , gap_MH, color = "tab:green")
# # plt.plot(x,exp.(model(x,paramMH)),linestyle = "dashed", label = label_scatterMH,color = "tab:green")



# # plt.title(L"MIS Gap Scaling Comparison HP = 2 deg = 1 $\beta = $"*string(beta))
# plt.ylabel(L"$\delta$")

# plt.xlabel(L"$N$")
# plt.yscale("log")
# # plt.xscale("log")
# plt.grid(which="major", zorder=-1.0)
# plt.legend()
# # name = "Figures/MIS/Classical/scalingBeta"*string(beta)*"HP2deg1.png"
# # name = "Figures/MIS/qHMC/scalingBeta"*string(beta)*"HP2deg1.png"
# # plt.savefig(name)
# plt.show()

