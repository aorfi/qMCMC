using Graphs
using Bloqade
using SparseArrays

function MIS_ham(atoms, ν)
    graph = BloqadeMIS.unit_disk_graph(atoms, 7.5)
    N = length(atoms)
    dim = 2^N
    H = zeros(dim,dim)
    for confi in (0:dim-1)
        Diagonal = 0
        for SpinIndex in (0:N-1)
            Si = 2*((confi>>( SpinIndex))&1)-1
            ni = 0.5*(1+Si)
            Diagonal += -ni
            neigh = neighbors(graph, SpinIndex+1) .- 1 
            for i in neigh
                Si_neigh = 2*((confi>>(i))&1)-1
                ni_neigh = 0.5*(1+Si_neigh)
                # This gets double counted so divide by 2
                Diagonal += ν*0.5*ni*ni_neigh 
            end
        end
        H[confi+1,confi+1] += Diagonal
    end
    return H |> sparse
end 