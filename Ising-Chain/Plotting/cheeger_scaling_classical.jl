using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2
using LsqFit
cd(pwd, "..")

N_values = (4:11) 
beta = 5

gap_glab = load_object("Data/Ising-Chain/Classical/GlaubScaling")
gap_glab_loc = load_object("Data/Ising-Chain/Classical/GlaubLocScaling")
gap_glab_locPBC = load_object("Data/Ising-Chain/Classical/GlaubLocScalingPBC")
gap_MH = load_object("Data/Ising-Chain/Classical/MHUniformScaling")
gap_MH_loc = load_object("Data/Ising-Chain/Classical/MHLocScaling")
gap_MH_locPBC = load_object("Data/Ising-Chain/Classical/MHLocScalingPBC")
# gap_qHMC = load_object("Data/Ising/qHMC/BestScaling")
# gap_rand = load_object("Data/Ising/qHMC/RandScaling")

phi_MH = load_object("Data/Ising-Chain/Classical/CheegerMHUniformScaling")
phi_MH_loc = load_object("Data/Ising-Chain/Classical/CheegerMHLocScaling")
phi_MH_locPBC = load_object("Data/Ising-Chain/Classical/CheegerMHLocScalingPBC")


@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]

fitG = curve_fit(model,N_values,log.(gap_glab),p0)
paramG = fitG.param
scaleG = -round(paramG[1]/log(2), digits=4)

fitGl = curve_fit(model,log.(N_values),log.(gap_glab_loc),p0)
paramGl = fitGl.param
scaleGl = round(paramGl[1], digits=4)

fitGlPBC = curve_fit(model,log.(N_values),log.(gap_glab_locPBC),p0)
paramGlPBC = fitGlPBC.param
scaleGlPBC = round(paramGlPBC[1], digits=4)

fitMH = curve_fit(model,N_values,log.(gap_MH),p0)
paramMH = fitMH.param
scaleMH = -round(paramMH[1]/log(2), digits=4)

fitMHl = curve_fit(model,log.(N_values),log.(gap_MH_loc),p0)
paramMHl = fitMHl.param
scaleMHl = round(paramMHl[1], digits=4)

fitMHlPBC = curve_fit(model,log.(N_values),log.(gap_MH_locPBC),p0)
paramMHlPBC = fitMHlPBC.param
scaleMHlPBC = round(paramMHlPBC[1], digits=4)

x = range(4,15, length= 1000)

# label_scatterMH = L"Uniform MH $2^{-kN}$ k= "*string(scaleMH)
# plt.scatter(N_values, gap_MH, color = "tab:green")
# plt.plot(x,exp.(model(x,paramMH)),linestyle = "dashed", label = label_scatterMH,color = "tab:green")
# plt.scatter(N_values ,  2*phi_MH, color = "black")
# plt.scatter(N_values ,  0.5*phi_MH.^2, color = "black")


label_scatterMHl = L"MH Local OBC $N^{b}$ b= "*string(scaleMHl)
plt.scatter(N_values ,  gap_MH_loc, color = "tab:orange")
plt.plot(x,exp.(model(log.(x),paramMHl)),linestyle = "dashed", label = label_scatterMHl,color = "tab:orange")
plt.scatter(N_values ,  2*phi_MH_loc, color = "black")
plt.scatter(N_values ,  0.5*phi_MH_loc.^2, color = "black")

function OBC(n,b)
    return (2*(exp(-2*b)-exp(-4*b))/n+exp(-4*b))
end
function OBC_scale(n,b)
    return 2exp(b*(n-1))/(2*(2*cosh(b))^(n-1)-2*exp(b*(n-1)))
end
plt.plot(x, [2*OBC(i,beta) for i in x], color = "black",linestyle = "dashed")
plt.plot(x, [0.5*OBC(i,beta)^2 for i in x], color = "black",linestyle = "dashed")
# plt.plot(x, [2*OBC(i,beta)*OBC_scale(i,beta) for i in x], color = "blue",linestyle = "dashed")


# label_scatterMHlPBC = L"MH Local PBC $N^{b}$ b= "*string(scaleMHlPBC)
# plt.scatter(N_values ,  gap_MH_locPBC, color = "tab:red")
# plt.plot(x,exp.(model(log.(x),paramMHlPBC)),linestyle = "dashed", label = label_scatterMHlPBC,color = "tab:red")
# plt.scatter(N_values ,  2*phi_MH_locPBC, color = "black")
# # plt.scatter(N_values ,  0.5*phi_MH_locPBC.^2, color = "black")
# function PBC(n,b)
#     return (2*(exp(-b*(-n+4)))/(((2*sinh(b))^n+(2*cosh(b))^n)-2*exp(n*b)))
# end
# plt.plot(x, [2*PBC(i,beta) for i in x], color = "black",linestyle = "dashed")
# plt.plot(x, [exp(-4*beta) for i in x], color = "black",linestyle = "dashed")

# label_scatterG = L"Uniform Glauber $2^{-kN}$ k= "*string(scaleG)
# plt.scatter(N_values , gap_glab, color = "tab:blue")
# plt.plot(x,exp.(model(x,paramG)),linestyle = "dashed", label = label_scatterG,color = "tab:blue")


# label_scatterGl = L"Glaubler Local OBC $N^{b}$ b= "*string(scaleGl)
# plt.scatter(N_values ,  gap_glab_loc, color = "tab:blue")
# plt.plot(x,exp.(model(log.(x),paramGl)),linestyle = "dashed", label = label_scatterGl,color = "tab:blue")

# label_scatterGlPBC = L"Glaubler Local PBC $N^{b}$ b= "*string(scaleGlPBC)
# plt.scatter(N_values ,  gap_glab_locPBC, color = "tab:green")
# plt.plot(x,exp.(model(log.(x),paramGlPBC)),linestyle = "dashed", label = label_scatterGlPBC,color = "tab:green")




plt.title( L"Ising Chain Local Scaling $\beta = $"*string(beta))
plt.ylabel(L"$\delta$")

plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid(which="major", zorder=-1.0)
plt.legend()
# name = "Figures/Ising-Chain/Classical/localScalingBeta"*string(beta)*".png"
# plt.savefig(name)
plt.show()


