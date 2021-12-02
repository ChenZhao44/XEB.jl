function max_space_ids(ec::NestedEinsum{OMEinsum.DynamicEinCode{Int}}, ids_dic)
    OMEinsum.isleaf(ec) && return Int[]
    max_ids = Int[]
    max_sc = -Inf
    for arg in ec.args
        new_max_ids = max_space_ids(arg, ids_dic)
        new_sc = space_complexity(new_max_ids, ids_dic)
        if new_sc > max_sc
            max_ids = new_max_ids
            max_sc = new_sc
        end
    end
    for ids_x in ec.eins.ixs
        if space_complexity(ids_x, ids_dic) > max_sc
            max_ids = ids_x
            max_sc = space_complexity(ids_x, ids_dic)
        end
    end
    if space_complexity(ec.eins.iy, ids_dic) > max_sc
        max_ids = ec.eins.iy
        max_sc = space_complexity(ec.eins.iy, ids_dic)
    end
    return copy(max_ids)
end
space_complexity(ids::Vector{Int}, ids_dic) = isempty(ids) ? -Inf : sum(log2(ids_dic[i]) for i in ids)
