using Graphs
using Bloqade
using Random
using LinearAlgebra
using SparseArrays



function mixing_rydberg_SA(H, atoms, beta, Rb)
    graph = BloqadeMIS.unit_disk_graph(atoms, Rb)
    N = length(atoms)
    dim = 2^N
    M = zeros(dim,dim)
    for row in (0:dim-1)
        E_sp = H[row+1,row+1]
        p_sum = 0 
        for SpinIndex in (0:N-1)
            # check if spin is up (in the set) or free
            in_set = mod(row,2^(SpinIndex+1))>= 2^(SpinIndex)
            neigh = neighbors(graph, SpinIndex+1) .- 1
            is_free = reduce(*,[!(mod(row,2^(i+1))>= 2^(i)) for i in neigh]) 
            if in_set
                num_neigh = length(neigh)
                # get binary of flip
                row_flip = row ⊻ (Int(2)^(SpinIndex))
                for i in neigh
                    # check if neihbour is also in the set 
                    neigh_in_set = mod(row,2^(i+1))>= 2^(i)
                    if !neigh_in_set
                        bit = Int(2)^(i)
                        col= row_flip ⊻ bit
                        E_s = H[col+1,col+1]
                        M[row+1,col+1] += 1/(N*8)*min(1,real(exp(beta*(E_sp-E_s))))
                        p_sum += 1/(N*8)*min(1,real(exp(beta*(E_sp-E_s))))
                    end
                end
                # flip it with the remaining probability
                remain_prob = 8-num_neigh
                bit = Int(2)^SpinIndex
                col = row ⊻ bit
                E_s = H[col+1,col+1]
                M[row+1,col+1] += remain_prob/(N*8)*min(1,real(exp(beta*(E_sp-E_s))))
                p_sum += remain_prob/(N*8)*min(1,real(exp(beta*(E_sp-E_s))))
            elseif is_free
                # flip the spin 
                bit = Int(2)^SpinIndex
                col = row ⊻ bit
                E_s = H[col+1,col+1]
                M[row+1,col+1] += 1/(N)*min(1,real(exp(beta*(E_sp-E_s))))
                p_sum += 1/(N)*min(1,real(exp(beta*(E_sp-E_s))))
            end
        end
        M[row+1,row+1] = 1-p_sum
    end
    return M |> sparse
end

function mixing_MIS_SA(H, ϵ, atoms, beta,Rb)
    graph = BloqadeMIS.unit_disk_graph(atoms, Rb)
    N = length(atoms)
    dim = 2^N
    M = zeros(dim,dim)
    for row in (0:dim-1)
        E_sp = H[row+1,row+1]
        p_sum = 0 
        for SpinIndex in (0:N-1)
            # check if spin is up (in the set)
            in_set = mod(row,2^(SpinIndex+1))>= 2^(SpinIndex)
            # get neighbours
            neigh = neighbors(graph, SpinIndex+1) .- 1 
            if in_set
                # remove with prob ϵ
                bit = Int(2)^SpinIndex
                col = row ⊻ bit 
                E_s = H[col+1,col+1]
                M[row+1,col+1] += (ϵ/N)*min(1,real(exp(beta*(E_sp-E_s))))
                p_sum += (ϵ/N)*min(1,real(exp(beta*(E_sp-E_s))))

                # spin exchange with neigh with prob (1-ϵ)/8
                # get binary of flip
                row_flip = row ⊻ (Int(2)^(SpinIndex))
                for i in neigh
                    # check if neihbour is also in the set 
                    neig_in_set = mod(row,2^(i+1))>= 2^(i)
                    if !neig_in_set
                        bit = Int(2)^(i)
                        col= row_flip ⊻ bit
                        E_s = H[col+1,col+1]
                        m = ((1-ϵ)/(8*N))*min(1,real(exp(beta*(E_sp-E_s))))
                        M[row+1,col+1] += m 
                        p_sum += m
                    end
                end
                # No update with the remaining probability is done in the p_sum
            else
                # flip the spin 
                bit = Int(2)^SpinIndex
                col = row ⊻ bit
                E_s = H[col+1,col+1]
                M[row+1,col+1] += (1/N)*min(1,real(exp(beta*(E_sp-E_s))))
                p_sum += (1/N)*min(1,real(exp(beta*(E_sp-E_s))))
            end
        end
        M[row+1,row+1] = 1-p_sum
    end
    return M |> sparse
end



