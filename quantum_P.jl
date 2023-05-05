
function mixing_ham(N) 
    dim = (2)^N
    Hm = zeros(dim,dim)
    for config in (0:dim-1)
        for spin_index in (0:N-1)
            binary_spin_index = Int(2)^(spin_index)
            config_flip = config  ⊻ binary_spin_index
            # Sigma_x term
            Hm[config_flip+1,config+1] += 1
        end
    end
    return Hm |> sparse
end

function mixing_qHMC(N, Hc,κ,η,beta, Hm)
    ratio = norm(real.(Matrix(Hm)),2)/norm(real.(Matrix(Hc)), 2)
    H = κ*ratio*Hc+η*Hm
    dim = 2^N
    M = zeros(dim,dim)
    U = exp(-1im*Matrix(H))
    prob = (abs.(U)).^2
    @inbounds for row in (0:dim-1)
        E_sp = Hc[row+1,row+1]
        diag = 0 
        for col in (0:dim-1)
            if row != col
                E_s = Hc[col+1,col+1]
                # MH
                m = prob[row+1,col+1]*min(1,real(exp(beta*(E_sp-E_s))))
                # Glaub
                # m = prob[row+1,col+1]*(1/(1+exp(beta*(E_s-E_sp))))
                M[row+1,col+1] = m
                diag += m
            end
        M[row+1,row+1] = 1-diag
        end
    end
    return M |> sparse
end
