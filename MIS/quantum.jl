using Bloqade
using Graphs

function mixing_exchange(atoms, Rb) 
    graph = BloqadeMIS.unit_disk_graph(atoms, Rb)
    N = length(atoms)
    dim = 2^N
    Hm = zeros(dim,dim)
    for row in (0:dim-1)
        for SpinIndex in (0:N-1)
            in_set = mod(row,2^(SpinIndex+1))>= 2^(SpinIndex)
            neigh = neighbors(graph, SpinIndex+1) .- 1
            if in_set
                row_flip = row ⊻ (Int(2)^(SpinIndex))
                for i in neigh
                    # check if neihbour is also in the set 
                    in_set = mod(row,2^(i+1))>= 2^(i)
                    if !in_set
                        bit = Int(2)^(i)
                        col= row_flip ⊻ bit
                        Hm[row+1,col+1] += 1
                    end
                end
            end
        end
    end
    return Hm
end


function mixing_rydberg_SA_qHMC(Hc, atoms,α,η, beta, Hm, Rb)
    graph = BloqadeMIS.unit_disk_graph(atoms, Rb)
    N = length(atoms)
    dim = 2^N
    M = zeros(dim,dim)
    ratio = norm(real.(Matrix(Hm)),2)/norm(real.(Matrix(Hc)), 2)
    H = α*ratio*Hc+η*Hm
    U = exp(-1im*Matrix(H))
    prob = (abs.(U)).^2
    for row in (0:dim-1)
        E_sp = Hc[row+1,row+1]
        p_sum = 0 
        for SpinIndex in (0:N-1)
            # check if spin is up (in the set) or free
            in_set = mod(row,2^(SpinIndex+1))>= 2^(SpinIndex)
            neigh = neighbors(graph, SpinIndex+1) .- 1
            is_free = reduce(*,[!(mod(row,2^(i+1))>= 2^(i)) for i in neigh]) 
            if in_set
                num_neigh = length(neigh)
                # get binary of flip
                for col in (0:dim-1)
                    if row != col
                        E_s = Hc[col+1,col+1]
                        m = prob[row+1,col+1]*(num_neigh/(8*N))*min(1,real(exp(beta*(E_sp-E_s))))
                        M[row+1,col+1] += m
                        p_sum += m 
                    end
                end
                # flip it with the remaining probability
                remain_prob = 8-num_neigh
                bit = Int(2)^SpinIndex
                col = row ⊻ bit
                E_s = Hc[col+1,col+1]
                M[row+1,col+1] += remain_prob/(N*8)*min(1,real(exp(beta*(E_sp-E_s))))
                p_sum += remain_prob/(N*8)*min(1,real(exp(beta*(E_sp-E_s))))
            elseif is_free
                # flip the spin 
                bit = Int(2)^SpinIndex
                col = row ⊻ bit
                E_s = Hc[col+1,col+1]
                M[row+1,col+1] += 1/(N)*min(1,real(exp(beta*(E_sp-E_s))))
                p_sum += 1/(N)*min(1,real(exp(beta*(E_sp-E_s))))
            end
        end
        M[row+1,row+1] = 0.9999999999999999-p_sum
    end
    return M # |> sparse
end

function mixing_MIS_SA_qHMC(Hc, ϵ, atoms,α,η, beta, Hm, Rb)
    graph = BloqadeMIS.unit_disk_graph(atoms, Rb)
    N = length(atoms)
    dim = 2^N
    M = zeros(dim,dim)
    ratio = norm(real.(Matrix(Hm)),2)/norm(real.(Matrix(Hc)), 2)
    H = α*ratio*Hc+η*Hm
    U = exp(-1im*Matrix(H))
    prob = (abs.(U)).^2
    for row in (0:dim-1)
        E_sp = Hc[row+1,row+1]
        p_sum = 0 
        for SpinIndex in (0:N-1)
            # check if spin is up (in the set)
            in_set = mod(row,2^(SpinIndex+1))>= 2^(SpinIndex)
            if in_set
                # remove with prob ϵ
                bit = Int(2)^SpinIndex
                col = row ⊻ bit 
                E_s = Hc[col+1,col+1]
                M[row+1,col+1] += (ϵ/N)*min(1,real(exp(beta*(E_sp-E_s))))
                p_sum += (ϵ/N)*min(1,real(exp(beta*(E_sp-E_s))))

                # evolve with prob (1-ϵ)
                for col in (0:dim-1)
                    if row != col
                        E_s = Hc[col+1,col+1]
                        m = prob[row+1,col+1]*(1-ϵ)/(N)*min(1,real(exp(beta*(E_sp-E_s))))
                        M[row+1,col+1] += m
                        p_sum += m 
                    end
                end
                # No update with the remaining probability is done in the p_sum
            else
                # flip the spin 
                bit = Int(2)^SpinIndex
                col = row ⊻ bit
                E_s = Hc[col+1,col+1]
                M[row+1,col+1] += (1/N)*min(1,real(exp(beta*(E_sp-E_s))))
                p_sum += (1/N)*min(1,real(exp(beta*(E_sp-E_s))))
            end
        end
        M[row+1,row+1] = 1-p_sum
    end
    return M |> sparse
end