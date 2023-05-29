using Bloqade
using PythonCall
using Arpack
using Random
using GenericTensorNetworks
using JLD2
plt = pyimport("matplotlib.pyplot")
# BloqadeLattices.DEFAULT_BACKGROUND_COLOR[] = "#FFFFFF"


function hardness_parameter(atoms,Rb)
    # Calculates hardness parameter
    graph = BloqadeMIS.unit_disk_graph(atoms, Rb)
    polynomial = GenericTensorNetworks.solve(IndependentSet(graph), GraphPolynomial())[]
    size_mis = length(polynomial)-1
    HP = polynomial[end-1]/(polynomial[end]*size_mis)
    return HP
end

function get_graph_mask(atoms, side_len, scale)
    grid = generate_sites(SquareLattice(), side_len, side_len; scale = scale)
    mask = zeros(side_len^2)
    [mask[i] = 1 for i in (1:side_len^2) if grid[i] in atoms] 
    mask = transpose(reshape(mask, (side_len, side_len)))
    return mask
end



a = 4.5
Rb = 7.5

# Generate 1000 configurations
# N_all = [4,5,6,7,8,10,11,13,15,18,20]
# L = [3,3,3,3,4,4,4,4,5,5,5]
# D_all = [0.6,0.4,0.3,0.2,0.5,0.4,0.3,0.2,0.4,0.3,0.2]
# for i in (1:11)
#     side_len = L[i]
#     D = D_all[i]
#     num = 1000
#     config_info = Vector{Any}(undef, num)
#     HP = zeros(num)
#     for i in (1:num)
#         println("working on ", i)
#         graph_index = i
#         atoms = generate_sites(SquareLattice(), side_len, side_len ; scale = a)|> random_dropout(D)
#         graph = BloqadeMIS.unit_disk_graph(atoms, Rb)
#         polynomial = GenericTensorNetworks.solve(IndependentSet(graph), GraphPolynomial())[]
#         size_mis = length(polynomial)-1
#         HP[i] = polynomial[end-1]/(polynomial[end]*size_mis)
#         graph_data = Dict("graph_index" => graph_index,"MIS_size" => size_mis,
#                             "degeneracy" => polynomial[end], "side_length" => side_len,
#                             "graph_mask" => get_graph_mask(atoms, side_len, a), "number_of_nodes" => length(atoms),
#                             "HP" => HP[i]) 
#         config_info[i] = graph_data
#         println("N = ",length(atoms))
#     end
#     save_object("Data/MIS/Configurations/L"*string(side_len)*"D"*string(D)*"config", config_info)
# end


# Get a configurations with deg = 1 and HP =2 for each N
N_all = [4,5,6,7,8,10,11,13,15,18,20]
L = [3,3,3,3,4,4,4,4,5,5,5]
D_all = [0.6,0.4,0.3,0.2,0.5,0.4,0.3,0.2,0.4,0.3,0.2]
num_N = length(N_all)
config_sorted = Vector{Any}(undef, num_N)
counter = 1
for i in (1:num_N)
    side_len = L[i]
    D = D_all[i]
    configs = load_object("Data/MIS/Configurations/L"*string(side_len)*"D"*string(D)*"config")
    num = length(configs)
    deg = [configs[i]["degeneracy"] for i in (1:num)]
    HP = [configs[i]["HP"] for i in (1:num)]
    for j in (1:num)
        if deg[j] == 1 && HP[j] == 2
            config_sorted[counter] = configs[i]
            println("N = ", N_all[i])
            global counter += 1
            break 
        end
    end
end
save_object("Data/MIS/Configurations/Deg1HP2configs", config_sorted)


