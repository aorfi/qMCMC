using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2
using LsqFit
cd(pwd, "..")

N_values = (4:11) 
beta = 4

# gap_glab = load_object("Data/Ising/Classical/GlaubScaling")
# gap_glab_loc = load_object("Data/Ising/Classical/GlaubLocScaling")
gap_MH = load_object("Data/Ising-Chain/Classical/MHUniformScaling")
gap_MH_loc = load_object("Data/Ising-Chain/Classical/MHLocScaling")
gap_MH_locPBC = load_object("Data/Ising-Chain/Classical/MHLocScalingPBC")
# gap_qHMC = load_object("Data/Ising/qHMC/BestScaling")
# gap_rand = load_object("Data/Ising/qHMC/RandScaling")

@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]

# fitq = curve_fit(model,N_values[1:end-1],log.(gap_qHMC),p0)
# paramq = fitq.param
# scaleq = -round(paramq[1]/log(2), digits=4)

# fitRand = curve_fit(model,N_values,log.(gap_rand),p0)
# paramRand = fitRand.param
# scaleRand = -round(paramRand[1]/log(2), digits=4)

# fitG = curve_fit(model,N_values,log.(gap_glab),p0)
# paramG = fitG.param
# scaleG = -round(paramG[1]/log(2), digits=4)

# fitGl = curve_fit(model,log.(N_values),log.(gap_glab_loc),p0)
# paramGl = fitGl.param
# scaleGl = round(paramGl[1], digits=4)

fitMH = curve_fit(model,N_values,log.(gap_MH),p0)
paramMH = fitMH.param
scaleMH = -round(paramMH[1]/log(2), digits=4)

fitMHl = curve_fit(model,log.(N_values),log.(gap_MH_loc),p0)
paramMHl = fitMHl.param
println(paramMHl)
scaleMHl = round(paramMHl[1], digits=4)

x = range(4,20, length= 1000)

label_scatterMH = L"Uniform MH $2^{-kN}$ k= "*string(scaleMH)
plt.scatter(N_values, gap_MH, color = "tab:green")
plt.plot(x,exp.(model(x,paramMH)),linestyle = "dashed", label = label_scatterMH,color = "tab:green")

label_scatterMHl = L"MH Local $N^{b}$ b= "*string(scaleMHl)
plt.scatter(N_values ,  gap_MH_loc, color = "tab:red")
plt.plot(x,exp.(model(log.(x),paramMHl)),linestyle = "dashed", label = label_scatterMHl,color = "tab:red")
plt.scatter(N_values ,  gap_MH_locPBC, color = "tab:blue")


# label_scatterG = L"Uniform Glauber $2^{-kN}$ k= "*string(scaleG)
# plt.scatter(N_values , gap_glab, color = "tab:blue")
# plt.plot(x,exp.(model(x,paramG)),linestyle = "dashed", label = label_scatterG,color = "tab:blue")


# label_scatterGl = L"Glaubler Local $N^{b}$ b= "*string(scaleGl)
# plt.scatter(N_values ,  gap_glab_loc, color = "tab:orange")
# plt.plot(x,exp.(model(log.(x),paramGl)),linestyle = "dashed", label = label_scatterGl,color = "tab:orange")


# label_scatterq  = L"qHMC Best $2^{-kN}$ k= "*string(scaleq)
# plt.scatter(N_values[1:end-1] ,  gap_qHMC, color = "tab:purple")
# plt.plot(x,exp.(model(x,paramq)),linestyle = "dashed", label = label_scatterq,color = "tab:purple")

# label_scatterRand  = L"qHMC Random $2^{-kN}$ k= "*string(scaleRand)
# plt.scatter(N_values ,  gap_rand, color = "tab:brown")
# plt.plot(x,exp.(model(x,paramRand)),linestyle = "dashed", label = label_scatterRand,color = "tab:brown")

plt.title( L"1D Ising OBC Scaling $\beta = $"*string(beta))
plt.ylabel(L"$\delta$")

plt.xlabel(L"$N$")
# plt.yscale("log")
# plt.xscale("log")
plt.grid(which="major", zorder=-1.0)
plt.legend()
# name = "Figures/Ising/qHMC/scalingBeta"*string(beta)*".png"
# plt.savefig(name)
plt.show()


