using LinearAlgebra
using Arpack
using JLD2
using SparseArrays

function ising_ham(N::Integer, couplings, h)
    dim = 2^N
    H = zeros(dim,dim)
    for config in (0:dim-1)
        Diagonal = 0
        for spin_index in (0:N-1)
            # Boundary condition
            if spin_index == N-1
                Si = 2*((config>>spin_index)&1)-1
                Diagonal += h*Si 
                Si_next = 2*((config>>(0))&1)-1
                Diagonal += -couplings[N]*Si*Si_next
                break
            end
            Si = 2*((config>>spin_index)&1)-1
            Diagonal += h*Si 
            Si_next = 2*((config>>(spin_index+1))&1)-1
            Diagonal += -couplings[N-(spin_index+1)]*Si*Si_next
        end
        H[config+1,config+1] += Diagonal
    end
    return H |> sparse
end 

function ising_energy(N::Integer, couplings, h, config)
    config in (0:2^N-1) || error("invalid configuration")
    eng = 0
    for spin_index in (0:N-1)
        if spin_index == N-1
            Si = 2*((config>>spin_index)&1)-1
            eng += h*Si
            Si_next = 2*((config>>(0))&1)-1
            eng += -couplings[N]*Si*Si_next
            break
        end
        Si = 2*((config>>spin_index)&1)-1
        eng += h*Si
        Si_next = 2*((config>>(spin_index+1))&1)-1
        eng += -couplings[N-(spin_index+1)]*Si*Si_next
    end
    return eng
end