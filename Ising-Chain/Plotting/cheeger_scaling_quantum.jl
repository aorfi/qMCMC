using GLM
using DataFrames
using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2

cd(pwd, "..")

N_values = (50:4:150)
beta = 5
phi_av = 2*load_object("Data/Ising-Chain/qMCMC/Cheeger/AverageBottleneckScaling")[:,1]
std_av = load_object("Data/Ising-Chain/qMCMC/Cheeger/AverageBottleneckScaling")[:,2]
# gap_glauber = load_object("Data/Ising-Chain/Classical/GlaubLocScalingAnalytic")
phi_lh = 2*load_object("Data/Ising-Chain/qMCMC/Cheeger/LowerHalfBottleneckScaling")[:,1]
std_lh = load_object("Data/Ising-Chain/qMCMC/Cheeger/LowerHalfBottleneckScaling")[:,2]
phi = 2*load_object("Data/Ising-Chain/qMCMC/Cheeger/BottleneckScalingEta40Kappa20")[:,1]

df = DataFrame(x = log.(N_values), y = log.(phi_av ))
y_values = log.(phi_av )
y_err = zeros(length(y_values))
for (i,sigma) in pairs(std_av)
    y_err[i] = abs(1/(y_values[i]))*sigma
end
df.w = 1 ./ y_err
fitq  = glm(@formula(y ~ x), df,Normal(), wts = df.w)
# display(fitq)
paramq = coef(fitq)
scaleq = -round(paramq[2], digits=4)
errorq = stderror(fitq)[2]
println(errorq)

x = range(50,150, length= 1000)

label_scatterq  = L"Average $\Phi$  $N^{-b}$ b= "*string(scaleq)*"(2)"
plt.errorbar(N_values ,  phi_av , yerr = std_av, fmt=".k",color = "tab:purple")
plt.plot(x,[exp(paramq[2]*log(i)+paramq[1]) for i in x],linestyle = "dashed", label = label_scatterq,color = "tab:purple")


x = N_values
y = log.(phi_lh )
fitq  = lm(@formula(y ~ x), (;x,y))
# display(fitq)
paramq = coef(fitq)
scaleq = -round(paramq[2]/log(2), digits=3)
errorq = stderror(fitq)[2]/log(2)
# println(errorq)


plt.errorbar(N_values ,  phi_lh , yerr = std_lh, fmt=".k",color = "tab:brown")
# plt.scatter(N_values, gap_glauber)
# plt.plot(x, [2^(-i+1) for i in x])
label_scatterq  = L"Average $\eta>\kappa$ $\Phi$ $2^{-kN}$ k= "*string(scaleq)*"(6)"
plt.plot(x,[exp(paramq[2]*i+paramq[1]) for i in x],linestyle = "dashed", label = label_scatterq,color = "tab:brown")

x = N_values
y = log.(phi)
fitq  = lm(@formula(y ~ x), (;x,y))
display(fitq)
paramq = coef(fitq)
scaleq = -round(paramq[2], digits=2)
scaleq2 = -round(paramq[2]/log(2), digits=2)
errorq = stderror(fitq)[2]/log(2)
println(errorq)


plt.scatter(N_values ,  phi ,color = "tab:red")
label_scatterq  = L" $\eta=40,\kappa=10$ $\Phi$ $2^{-kN}$ k= "*string(scaleq2)*L"(4), $e^{-\lambda N}$ $\lambda=$ "*string(scaleq)*"(3)"
plt.plot(x,[exp(paramq[2]*i+paramq[1]) for i in x],linestyle = "dashed", label = label_scatterq,color = "tab:red")

plt.ylabel(L"$\Phi$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid(which="major", zorder=-1.0)
# plt.ylim(10^(-16),1)
plt.legend()
# name = "Figures/Ising-Chain/scalingBeta"*string(beta)*".svg"
# plt.savefig(name)
plt.show()


