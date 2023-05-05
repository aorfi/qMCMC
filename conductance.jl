using IterTools

function min_conductance_search(Hc, P, N, beta)
    energies = [Hc[i,i] for i in (1:2^N)]
    steady_state = normalize([exp(-beta*e) for e in energies],1)
    possible_configs = (1:2^N)
    all_subsets = collect(subsets(possible_configs))[2:end]
    all_phi_A = Float64[]
    subset_small = Any[]
    for i in (1:length(all_subsets))
        set = all_subsets[i]
        pi_set = 0 
        for i in set
            pi_set += steady_state[i]
        end
        if pi_set <= 0.5
            Q_A = 0
            for i in set
                Q_i = 0
                for j in (1:2^N)
                    if !(j in set)
                        Q_i += steady_state[i]*P[i,j]
                    end
                end
                Q_A += Q_i
            end
            phi_A = Q_A/pi_set
            append!(all_phi_A, phi_A)
            append!(subset_small, [set])
        end
    end
    # println(all_phi_A)
    # println(subset_small)
    phi_min = findmin(all_phi_A)[1]
    min_A = subset_small[findmin(all_phi_A)[2]]
    return phi_min, min_A, all_phi_A,subset_small
end

function conductance(Hc, P, N, A,beta)
    energies = [Hc[i,i] for i in (1:2^N)]
    steady_state = normalize([exp(-beta*e) for e in energies],1)
    set = A
    pi_set = 0 
    Q_A = 0
    for i in set
        pi_set += steady_state[i]
        Q_i = 0
        for j in (1:2^N)
            if !(j in set)
                Q_i += steady_state[i]*P[i,j]
            end
        end
        Q_A += Q_i
    end
    phi_A = Q_A/pi_set
    return phi_A
end


