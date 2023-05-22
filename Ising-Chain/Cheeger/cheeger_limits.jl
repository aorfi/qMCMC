using QuadGK
using PythonCall
plt = pyimport("matplotlib.pyplot")

function λ_sum(N, κ,η)
    k_pos0 = (2*pi/N)*range(1,(N)/2) .- pi/N
    overlap = []
    for k in k_pos0
        ek = sqrt((η-κ*cos(k))^2+(κ*sin(k))^2)
        append!(overlap, log(1- η^2*sin(k)^2*sin(2*ek)^2/ek^2) )
    end
    λ = -sum(overlap)/N
    return λ
end

function λ_int(k, α,η)
    ek = sqrt((η-α*cos(k))^2+(α*sin(k))^2)
    return -log(1- η^2*sin(k)^2*sin(2*ek)^2/ek^2)/(2*π)
end 

function γ_int(k, α,η)
    ek = sqrt((η-α*cos(k))^2+(α*sin(k))^2)
    return η^2*sin(k)^2*sin(2*ek)^2/((ek^2-η^2*sin(k)^2*sin(2*ek)^2)*(2*π))
end 

# function κ_int(k, α,η) # rename this!
#     ek = sqrt((η-α*cos(k))^2+(α*sin(k))^2)
#     return -log(cos(2*ek)- 1im*(η*cos(k)-α)*sin(2*ek)/ek)/(2*π)
# end 


# κ = 10
# num_values = 100
# eta_values = range(0.00000001,30, length=num_values)
# integral_λ = zeros(num_values)
# integral_γ_full = zeros(num_values)
# integral_γ = zeros(num_values)
# for eta_i in (1:num_values)
#     η = eta_values[eta_i]
#     integral_λ[eta_i], error = quadgk(k -> λ_int(k,κ,η), 0, π)
#     integral_γ_full[eta_i], error = quadgk(k -> γ_int(k,κ,η), 0, π)
#     k_values = range(0.0000000001,π, length=num_values)
#     γ = []
#     for k in k_values
#         ek = sqrt((η-κ*cos(k))^2+(κ*sin(k))^2) 
#         append!(γ, η^2*sin(k)^2*sin(2*ek)^2/((ek^2-η^2*sin(k)^2*sin(2*ek)^2)*(2π)))
#     end
#     integral_γ[eta_i] = sum(γ.*(k_values[2]-k_values[1]))
# end
# # plt.plot(eta_values,integral_λ)
# plt.plot(eta_values,integral_γ)
# plt.plot(eta_values,integral_γ_full)
# plt.ylim(0,20)
# plt.show()

κ = 5
l = 20
η = l*κ
k_values = range(0.0000000001,π, length=num_values)
test = []
test1 = []
test2 = []
for k in k_values
    ek = sqrt((η-κ*cos(k))^2+(κ*sin(k))^2) 
    append!(test,1- η^2*sin(k)^2*sin(2*ek)^2/ek^2)
    append!(test1,1- η^2*sin(k)^2/ek^2)
    append!(test2,1-sin(k)^2*sin(2*ek)^2)
end
int = -log.(test)/(2π)
result = sum(int.*(k_values[2]-k_values[1]))
int1 = -log.(test1)/(2π)
result1 = sum(int1.*(k_values[2]-k_values[1]))
int2 = -log.(test2)/(2π)
result2 = sum(int2.*(k_values[2]-k_values[1]))
println(result)
println(result1)
println(result2)
println(result/result1)



κ = 10
num_values = 100
eta_values = range(0.00000001,50, length=num_values)
results = []
results1 = []
results2 = []
ratios = []
for eta_i in (1:num_values)
    η = eta_values[eta_i]
    k_values = range(0.0000000001,π, length=num_values)
    test = []
    test1 = []
    test2 = []
    for k in k_values
        ek = sqrt((η-κ*cos(k))^2+(κ*sin(k))^2) 
        append!(test,1- η^2*sin(k)^2*sin(2*ek)^2/ek^2)
        append!(test1,1- η^2*sin(k)^2/ek^2)
        append!(test2,1-sin(k)^2*sin(2*ek)^2)
    end
    int = -log.(test)/(2π)
    result = sum(int.*(k_values[2]-k_values[1]))
    int1 = -log.(test1)/(2π)
    result1 = sum(int1.*(k_values[2]-k_values[1]))
    int2 = -log.(test2)/(2π)
    result2 = sum(int2.*(k_values[2]-k_values[1]))
    println(result/result1)
    append!(results,result)
    append!(results1,result1)
    append!(results2,result2)
    append!(ratios,result/result1) 
end
# plt.plot(eta_values,ratios)
plt.plot(eta_values,results)
plt.plot(eta_values,results1)
plt.plot(eta_values,results2)
plt.legend()
plt.show()




# N=50
# num_values = 100
# kappa_values = range(0.00000001,30, length=num_values)
# eta_values = range(0.00000001,30, length=num_values)
# integral_λ = zeros(num_values,num_values)
# integral_γ = zeros(num_values,num_values)
# for kappa_i in (1:num_values)
#     κ = kappa_values[kappa_i]
#     # println("κ ", κ)
#     for eta_i in (1:num_values)
#         η = eta_values[eta_i]
#         integral_λ[kappa_i,eta_i], error= quadgk(k -> λ_int(k,κ,η), 0, π)
#         integral_γ[kappa_i,eta_i], error= quadgk(k -> γ_int(k,κ,η), 0, π)
#         # integral_λ[kappa_i,eta_i] = λ_sum(N, κ,η)
#     end
# end
# plt.imshow(integral_λ, origin="lower",extent = [0, 50 , 0, 50],aspect="auto")
# # plt.imshow(integral_γ, origin="lower",extent = [0, 30 , 0, 30],aspect="auto")
# # plt.imshow(exp.(-N*integral_λ).*integral_γ*2/(N-1), origin="lower",extent = [0, 30 , 0, 30],aspect="auto")
# bar = plt.colorbar()
# # plt.clim(0, 1) 
# plt.xlabel(L"$\eta$")
# plt.ylabel(L"$\kappa$")
# bar.set_label(L"$λ$")
# plt.show()


