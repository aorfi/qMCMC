using GLM
using DataFrames
using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2

cd(pwd, "..")

N_values = (4:14) 
beta = 5


gap_MH = load_object("Data/Ising-Chain/Classical/MHUniformScaling")
gap_MH_loc = load_object("Data/Ising-Chain/Classical/MHLocScaling")
# gap_qHMC = load_object("Data/Ising-Chain/qMCMC/LargestGapScaling")
gap_av = load_object("Data/Ising-Chain/qMCMC/AverageGapScaling")[:,1]
std_av = load_object("Data/Ising-Chain/qMCMC/AverageGapScaling")[:,2]

x = N_values
y = log.(gap_MH) 
fitMH  = lm(@formula(y ~ x), (;x,y))
# display(fitMH)
paramMH = coef(fitMH)
scaleMH = -round(paramMH[2]/log(2), digits=2)

x = log.(N_values)
y = log.(gap_MH_loc) 
fitMHl  = lm(@formula(y ~ x), (;x,y))
# display(fitMHl)
paramMHl = coef(fitMHl)
scaleMHl = -round(paramMHl[2], digits=2)



df = DataFrame(x = N_values[1:end-6], y = log.(gap_av))
y_values = log.(gap_av)
y_err = zeros(length(y_values))
for (i,sigma) in pairs(std_av)
    y_err[i] = abs(1/(y[i]))*sigma
end
df.w = 1 ./ y_err
fitq  = glm(@formula(y ~ x), df,Normal(), wts = df.w)
display(fitq)
paramq = coef(fitq)
scaleq = -round(paramq[2]/log(2), digits=3)
errorq = stderror(fitq)[2]/log(2)
println(errorq)

x = range(4,15, length= 1000)

label_scatterMH = L"Uniform $2^{-kN}$ k= "*string(scaleMH)
plt.scatter(N_values, gap_MH, color = "tab:green")
plt.plot(x,[exp(paramMH[2]*i+paramMH[1]) for i in x],linestyle = "dashed", label = label_scatterMH,color = "tab:green")

label_scatterMHl = L"Local $N^{b}$ b= "*string(scaleMHl)
plt.scatter(N_values,  gap_MH_loc, color = "tab:red")
plt.plot(x,[exp(paramMHl[2]*log(i)+paramMHl[1]) for i in x],linestyle = "dashed", label = label_scatterMHl,color = "tab:red")

label_scatterq  = L"Average qMCMC $2^{-kN}$ k= "*string(scaleq)*"(3)"
plt.errorbar(N_values[1:end-6] ,  gap_av, yerr = std_av, fmt=".k",color = "tab:purple")
plt.plot(x,[exp(paramq[2]*i+paramq[1]) for i in x],linestyle = "dashed", label = label_scatterq,color = "tab:purple")


plt.ylabel(L"$\delta$")

plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid(which="major", zorder=-1.0)
plt.legend()
name = "Figures/Ising-Chain/scalingBeta"*string(beta)*".svg"
plt.savefig(name)
plt.show()


