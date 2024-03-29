function brick_layout(nbits::VT, ncycles::Integer) where {VT}
    edge_patterns = Dict{Edge{VT}, Symbol}()
    pattern_loop = [:A, :B]
    
    for i = 1:nbits-1
        if i % 2 == 0
            edge_patterns[Edge(i, i+1)] = :B
        else
            edge_patterns[Edge(i, i+1)] = :A
        end
    end
    qubit_locs = Dict{VT, Tuple{Float64, Float64}}()
    
    layout = RQCLayout{VT, Symbol, SCGate}(nbits, edge_patterns, pattern_loop, qubit_locs)
    if ncycles >= 1
        for i = 1:ncycles
            for v in vertices(layout)
                push!(gates(layout, v), Single)
                if has_partner(layout, v, i)
                    push!(gates(layout, v), FSim)
                else
                    push!(gates(layout, v), Id)
                end
            end
        end
    end
    for v in vertices(layout)
        push!(gates(layout, v), Single)
    end
    return layout
end

function google_layout_53(nbits::VT = 53, ncycles::Integer = 20) where {VT}
    edge_patterns = Dict{Edge{VT}, Symbol}()
    pattern_A = ((31,32),(21,22),(8,11),(24,29),(7,18),(1,5),(26,40),(15,25),
        (6,16),(44,53),(42,48),(46,51),(4,14),(13,27),(28,39),(2,3),(9,17),
        (23,30),(10,12),(19,20),(33,34),(41,47),(43,50),(45,49))
    pattern_B = ((32,37),(21,24),(18,26),(25,44),(22,35),(7,8),(5,15),(16,42),
        (1,4),(2,6),(12,51), (14,36),(3,13),(9,10),(20,41),(27,38),(17,28),
        (19,23),(34,43))
    pattern_C = ((32,52),(24,31),(26,29),(40,44),(22,37),(7,21),(15,18),(25,42),
        (11,35),(1,8),(5,6),(16,51),(3,4),(2,10),(12,41),(27,36),(13,17),(9,19),
        (20,43),(38,39),(28,30),(23,33),(34,45))
    pattern_D = ((21,32),(18,24),(25,26),(44,48),(8,22),(5,7),(15,16),(42,46),
        (4,11),(1,2),(6,12),(47,51),(13,14),(3,9),(10,20),(41,50),(27,28),
        (17,23),(19,34),(43,49))
    pattern_loop = [:A, :B, :C, :D, :C, :D, :A, :B]
    
    for e in pattern_A
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :A
    end
    for e in pattern_B
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :B
    end
    for e in pattern_C
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :C
    end
    for e in pattern_D
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :D
    end
    
    qubit_locs = Dict{VT, Tuple{Float64, Float64}}(
        53 => (0.0, 0), 48 => (2, 0), 46 => (4, 0), 47 => (6, 0), 
        50 => (8, 0), 49 => (10, 0), 44 => (1, 1), 42 => (3, 1), 
        51 => (5, 1), 41 => (7, 1), 43 => (9, 1), 45 => (11, 1),
        40 => (0, 2), 25 => (2, 2), 16 => (4, 2), 12 => (6, 2),
        20 => (8, 2), 34 => (10, 2), 26 => (1, 3), 15 => (3, 3),
        6 => (5, 3), 10 => (7, 3), 19 => (9, 3), 33 => (11, 3),
        29 => (0, 4), 18 => (2, 4), 5 => (4, 4), 2 => (6, 4),
        9 => (8, 4), 23 => (10, 4), 24 => (1, 5), 7 => (3, 5),
        1 => (5, 5), 3 => (7, 5), 17 => (9, 5), 30 => (11, 5),
        31 => (0, 6), 21 => (2, 6), 8 => (4, 6), 4 => (6, 6),
        13 => (8, 6), 28 => (10, 6), 32 => (1, 7), 22 => (3, 7),
        11 => (5, 7), 14 => (7, 7), 27 => (9, 7), 39 => (11, 7),
        52 => (0, 8), 37 => (2, 8), 35 => (4, 8), 36 => (8, 8),
        38 => (10, 8),
    )
    
    layout = RQCLayout{VT, Symbol, SCGate}(nbits, edge_patterns, pattern_loop, qubit_locs)
    if ncycles >= 1
        for i = 1:ncycles
            for v in vertices(layout)
                push!(gates(layout, v), Single)
                if has_partner(layout, v, i)
                    push!(gates(layout, v), FSim)
                else
                    push!(gates(layout, v), Id)
                end
            end
        end
    end
    for v in vertices(layout)
        push!(gates(layout, v), Single)
    end
    return layout
end

function square_layout(nbits::VT = 53, ncycles::Integer = 20, width = 10; start_with_zero = false) where {VT}
    edge_patterns = Dict{Edge{VT}, Symbol}()
    pattern_A = Tuple{Int,Int}[]
    pattern_B = Tuple{Int,Int}[]
    pattern_C = Tuple{Int,Int}[]
    pattern_D = Tuple{Int,Int}[]
    qubit_locs = Dict{VT, Tuple{Float64, Float64}}()
    if width % 2 == 0
        for i = 0:(nbits ÷ width)
            for j = 0:(width-1)
                if i == 0
                    j % 2 == 0 && (start_with_zero || j > 0) && push!(pattern_B, (j, j+1))
                    j % 2 == 1 && j < (width-1) && push!(pattern_D, (j, j+1))
                elseif i % 2 == 0
                    j % 2 == 0 && push!(pattern_B, (i*width-j-2, i*width+j))
                    j % 2 == 0 && j > 0 && push!(pattern_D, (i*width-j, i*width+j))
                    j % 2 == 1 && push!(pattern_B, (i*width-j, i*width+j))
                    j % 2 == 1 && j < (width-1) && push!(pattern_D, (i*width-j-2, i*width+j))
                elseif i % 2 == 1
                    j % 2 == 0 && push!(pattern_C, (i*width-j-2, i*width+j))
                    j % 2 == 0 && j > 0 && push!(pattern_A, (i*width-j, i*width+j))
                    j % 2 == 1 && push!(pattern_C, (i*width-j, i*width+j))
                    j % 2 == 1 && j < (width-1) && push!(pattern_A, (i*width-j-2, i*width+j))
                end
                i % 2 == 0 && j % 2 == 0 && (qubit_locs[i*width+j] = (j, i+0.5))
                i % 2 == 0 && j % 2 == 1 && (qubit_locs[i*width+j] = (j, -i-0.5))
                i % 2 == 1 && j % 2 == 0 && (qubit_locs[i*width+j] = (width-j-1, i+0.5))
                i % 2 == 1 && j % 2 == 1 && (qubit_locs[i*width+j] = (width-j-1, -i-0.5))
            end
        end
    else
        for i = 0:(nbits ÷ width)
            for j = 0:(width-1)
                if i == 0
                    j % 2 == 0 && j < (width-1) && (start_with_zero || j > 0) && push!(pattern_B, (j, j+1))
                    j % 2 == 1 && push!(pattern_D, (j, j+1))
                elseif i % 2 == 0
                    j % 2 == 0 && j < (width-1) && push!(pattern_B, (i*width-j-2, i*width+j))
                    j % 2 == 0 && j > 0 && push!(pattern_D, (i*width-j, i*width+j))
                    j % 2 == 1 && push!(pattern_B, (i*width-j, i*width+j))
                    j % 2 == 1 && j < (width-1) && push!(pattern_D, (i*width-j-2, i*width+j))
                elseif i % 2 == 1
                    j % 2 == 0 && j > 0 && push!(pattern_C, (i*width-j, i*width+j))
                    j % 2 == 0 && j < (width-1) && push!(pattern_A, (i*width-j-2, i*width+j))
                    j % 2 == 1 && push!(pattern_C, (i*width-j-2, i*width+j))
                    j % 2 == 1 && push!(pattern_A, (i*width-j, i*width+j))
                end
                i % 2 == 0 && j % 2 == 0 && (qubit_locs[i*width+j] = (j, i+0.5))
                i % 2 == 0 && j % 2 == 1 && (qubit_locs[i*width+j] = (j, -i-0.5))
                i % 2 == 1 && j % 2 == 0 && (qubit_locs[i*width+j] = (width-j-1, -i-0.5))
                i % 2 == 1 && j % 2 == 1 && (qubit_locs[i*width+j] = (width-j-1, i+0.5))
            end
        end
    end
    
    pattern_loop = [:A, :B, :C, :D, :C, :D, :A, :B]
    
    for e in pattern_A
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :A
    end
    for e in pattern_B
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :B
    end
    for e in pattern_C
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :C
    end
    for e in pattern_D
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :D
    end
    
    bits = start_with_zero ? (0:nbits-1) : (1:nbits)
    layout = RQCLayout{VT, Symbol, SCGate}(collect(bits), edge_patterns, pattern_loop, qubit_locs)
    if ncycles >= 1
        for i = 1:ncycles
            for v in vertices(layout)
                push!(gates(layout, v), Single)
                if has_partner(layout, v, i)
                    push!(gates(layout, v), FSim)
                else
                    push!(gates(layout, v), Id)
                end
            end
        end
    end
    for v in vertices(layout)
        push!(gates(layout, v), Single)
    end
    return layout
end

function ustc_layout_60(nbits = 60, ncycles = 24)
    vs = collect(0:65)
    vs_disable = [4, 5, 11, 17, 22, 65]
    vs_enable = setdiff(vs, vs_disable)[1:nbits]
    edge_patterns = Dict{Edge{Int}, Symbol}()
    pattern_A = ((0, 6), (1, 7), (2, 8), (9, 14), (12, 18), (13, 19), (15, 21), (20, 25), 
        (24, 30), (26, 32), (27, 33), (28, 34), (31, 36), (35, 40), (37, 43), (38, 44), 
        (39, 45), (41, 47), (48, 54), (49, 55), (50, 56), (52, 58), (53, 59), (57, 62))
    pattern_B = ((3, 9), (7, 12), (8, 13), (10, 15), (14, 20), (19, 24), (21, 26), (23, 28),
        (25, 31), (29, 35), (32, 37), (33, 38), (34, 39), (36, 42), (40, 46), (43, 48),
        (44, 49), (45, 50), (47, 52), (51, 57), (55, 60), (56, 61), (58, 63), (59, 64))
    pattern_C = ((1, 8), (6, 12), (7, 13), (9, 15), (10, 16), (14, 21), (18, 24), (19, 25), (20, 26), 
        (23, 29), (27, 34), (30, 36), (31, 37), (32, 38), (33, 39), (35, 41), (40, 47), (42, 48), 
        (43, 49), (44, 50), (45, 51), (46, 52), (54, 60), (55, 61), (56, 62), (57, 63), (58, 64))
    pattern_D = ((0, 7), (2, 9), (3, 10), (8, 14), (12, 19), (13, 20), (16, 23), (21, 27), 
        (24, 31), (25, 32), (26, 33), (28, 35), (34, 40), (36, 43), (37, 44), (38, 45), 
        (39, 46), (47, 53), (48, 55), (49, 56), (50, 57), (51, 58), (52, 59))
    pattern_loop = [:A, :B, :C, :D, :C, :D, :A, :B]
    
    for e in pattern_A
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :A
    end
    for e in pattern_B
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :B
    end
    for e in pattern_C
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :C
    end
    for e in pattern_D
        e[1] < e[2] || error("$e")
        edge_patterns[Edge(e)] = :D
    end
    
    layout = RQCLayout{Int, Symbol, SCGate}(vs_enable, edge_patterns, pattern_loop)
    if ncycles >= 1
        for i = 1:ncycles
            for v in vertices(layout)
                push!(gates(layout, v), Single)
                if has_partner(layout, v, i)
                    push!(gates(layout, v), FSim)
                else
                    push!(gates(layout, v), Id)
                end
            end
        end
    end
    for v in vertices(layout)
        push!(gates(layout, v), Single)
    end
    return layout
end