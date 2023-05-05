
function MH_uniform(N::Integer,H, beta)
    dim = (2)^N
    P = zeros(dim,dim)
    for row in (0:dim-1)
        E_sp = H[row+1,row+1]
        diag = 0 
        for col in (0:dim-1)
            if col != row
                E_s = H[col+1,col+1]
                m = (1/2^N)*min(1,real(exp(beta*(E_sp-E_s))))
                P[row+1,col+1] = m
                diag += m
            end
        P[row+1,row+1] = 1-diag
        end
    end
    return P
end

function MH_local(N::Integer,H, beta)
    dim = 2^N
    P = zeros(dim,dim)
    for config in (0:dim-1)
        E_sp = H[config+1,config+1]
        diag = 0
        for spin_index in (0:N-1)
            binary_spin_index = Int(2)^(spin_index)
            # finds config that differs by a spin flip
            config_flip = config ⊻ binary_spin_index
            E_s = H[config_flip+1,config_flip+1]
            m = (1/N)*min(1,real(exp(beta*(E_s-E_sp))))
            P[config_flip+1,config+1] = m
            diag += (1/N)*min(1,real(exp(-beta*(E_s-E_sp))))
        end 
        P[config+1,config+1] = 1-diag
    end
    return P |> sparse
end

function glab_uniform(N::Integer, H, beta)
    dim = 2^N
    P = zeros(dim,dim)
    for row in (0:dim-1)
        E_sp = H[row+1,row+1]
        diag = 0 
        for col in (0:dim-1)
            if col != row
                E_s = H[col+1,col+1]
                m = (1/2^N)*(1/(1+exp(beta*(E_s-E_sp))))
                P[row+1,col+1] = m
                diag += m
            end
        P[row+1,row+1] = 1-diag
        end
    end
    return P
end

function glab_local(N,H, beta)
    dim = 2^N
    P = zeros(dim,dim)
    for config in (0:dim-1)
        E_sp = H[config+1,config+1]
        diag = 0
        for spin_index in (0:N-1)
            binary_spin_index = Int(2)^(spin_index)
            # finds config that differs by a spin flip
            config_flip = config ⊻ binary_spin_index
            E_s = H[config_flip+1,config_flip+1]
            m = (1/N)*(1/(1+exp(-beta*(E_s-E_sp))))
            P[config_flip+1,config+1] = m
            diag += (1/N)*(1/(1+exp(beta*(E_s-E_sp))))
        end 
        P[config+1,config+1] = 1-diag
    end
    return P |> sparse
end
