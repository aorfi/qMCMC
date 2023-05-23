using LinearAlgebra
using Arpack
using SparseArrays
using JLD2

function SK_ham(N::Integer, couplings, h)
    # Fully connected
    dim = 2^N
    H = zeros(dim,dim)
    for config in (0:dim-1)
        Diagonal = 0
        for SpinIndex in (0:N-1)
            Si = 2*((config>>( SpinIndex))&1)-1
            Diagonal += h[SpinIndex+1]*Si 
            neigh = (SpinIndex+1:N-1)
            for i in neigh
                Si_neigh = 2*((config>>(i))&1)-1
                Diagonal += -couplings[sum(1:N-1)+1-SpinIndex-i]*Si*Si_neigh
            end
        end
        H[config+1,config+1] += Diagonal
    end
    return H |> sparse
end 

