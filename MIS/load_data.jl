using Bloqade


# Given a graph mask apply it to a Bloqade atoms object
function apply_graph_mask(sites::AtomList{D,T}, graph_mask) where {D,T}
    num_sites = size(sites)[1]
    mask_flat = reshape(transpose(graph_mask), num_sites)
    atoms = [i for i in 1:num_sites if mask_flat[i] == 1] 
    return sites[sort!(atoms)]
end

