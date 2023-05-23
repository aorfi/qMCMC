using GLM
using DataFrames
using PythonCall
plt = pyimport("matplotlib.pyplot")
using LaTeXStrings
using JLD2

cd(pwd, "..")

N_values = (4:10) 
beta = 5


gap_MH = load_object("Data/SK/Classical/MHScaling")[:,1]
gapMH_std = load_object("Data/SK/Classical/MHScaling")[:,2]
gap_MHl = load_object("Data/SK/Classical/MHLocScaling")[:,1]
gapMHl_std = load_object("Data/SK/Classical/MHLocScaling")[:,2]


x = range(4,15, length= 1000)


df = DataFrame(x = N_values, y = log.(gap_MH))
y_values = log.(gap_MH)
y_err = zeros(length(y_values))
for (i,sigma) in pairs(gapMH_std)
    y_err[i] = abs(1/(y_values[i]))*sigma
end
df.w = 1 ./ y_err
fit_MH  = glm(@formula(y ~ x), df,Normal(), wts = df.w)
display(fit_MH)
param_MH = coef(fit_MH)
scale_MH = -round(param_MH[2]/log(2), digits=3)
error_MH = stderror(fit_MH)[2]/log(2)
println(error_MH)

label_scatterMH = L"Uniform Proposal $2^{-kN}$ k= "*string(scale_MH)
plt.errorbar(N_values,  gap_MH, yerr = gapMH_std , fmt=".k",color = "tab:green")
plt.plot(x,[exp(param_MH[2]*i+param_MH[1]) for i in x],linestyle = "dashed", label = label_scatterMH,color = "tab:green")


# label_scatterMHl = L"Local Proposal $N^{-b}$ b= "*string(scaleMHl)
# plt.errorbar(N_values,  gap_MHl, yerr = gapMH_std , fmt=".k",color = "tab:blue")
# plt.plot(x,[exp(paramMHl[2]*log(i)+paramMHl[1]) for i in x],linestyle = "dashed", label = label_scatterMHl,color = "tab:blue")


plt.ylabel(L"$\delta$")

plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid(which="major", zorder=-1.0)
plt.legend()
name = "Figures/SK/scalingBeta"*string(beta)*".svg"
plt.savefig(name)
plt.show()


