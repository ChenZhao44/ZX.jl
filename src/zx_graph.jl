module EdgeType
    @enum EType SIM HAD
end

"""
    ZXGraph{T, P}

This is the type for representing the graph-like ZX-diagrams.
"""
struct ZXGraph{T<:Integer, P} <: AbstractZXDiagram{T, P}
    mg::Multigraph{T}
    ps::Dict{T, P}
    st::Dict{T, SpiderType.SType}
    et::Dict{Tuple{T, T}, EdgeType.EType}
    layout::ZXLayout{T}
    phase_ids::Dict{T,Tuple{T, Int}}
    scalar::Scalar{P}
    master::ZXDiagram{T, P}
    inputs::Vector{T}
    outputs::Vector{T}
end

function Base.copy(zxg::ZXGraph{T, P}) where {T, P}
    ZXGraph{T, P}(
        copy(zxg.mg), copy(zxg.ps),
        copy(zxg.st), copy(zxg.et), copy(zxg.layout),
        deepcopy(zxg.phase_ids), copy(zxg.scalar), 
        copy(zxg.master), copy(zxg.inputs), copy(zxg.outputs)
    )
end

"""
    ZXGraph(zxd::ZXDiagram)

Convert a ZX-diagram to graph-like ZX-diagram.

```jldoctest
julia> using ZXCalculus

julia> zxd = ZXDiagram(2); push_gate!(zxd, Val{:CNOT}(), 2, 1);

julia> zxg = ZXGraph(zxd)
ZX-graph with 6 vertices and 5 edges:
(S_1{input} <-> S_5{phase = 0//1⋅π})
(S_2{output} <-> S_5{phase = 0//1⋅π})
(S_3{input} <-> S_6{phase = 0//1⋅π})
(S_4{output} <-> S_6{phase = 0//1⋅π})
(S_5{phase = 0//1⋅π} <-> S_6{phase = 0//1⋅π})

```
"""
function ZXGraph(zxd::ZXDiagram{T, P}) where {T, P}
    zxd = copy(zxd)
    nzxd = copy(zxd)

    simplify!(Rule{:i1}(), nzxd)
    simplify!(Rule{:h}(), nzxd)
    simplify!(Rule{:i2}(), nzxd)
    match_f = match(Rule{:f}(), nzxd)
    while length(match_f) > 0
        for m in match_f
            vs = m.vertices
            if check_rule(Rule{:f}(), nzxd, vs)
                rewrite!(Rule{:f}(), nzxd, vs)
                v1, v2 = vs
                set_phase!(zxd, v1, phase(zxd, v1) + phase(zxd, v2))
                set_phase!(zxd, v2, zero(P))
            end
        end
        match_f = match(Rule{:f}(), nzxd)
    end

    vs = spiders(nzxd)
    vH = T[]
    vZ = T[]
    vB = T[]
    for v in vs
        if spider_type(nzxd, v) == SpiderType.H
            push!(vH, v)
        elseif spider_type(nzxd, v) == SpiderType.Z
            push!(vZ, v)
        else
            push!(vB, v)
        end
    end
    eH = [(neighbors(nzxd, v, count_mul = true)[1], neighbors(nzxd, v, count_mul = true)[2]) for v in vH]

    rem_spiders!(nzxd, vH)
    et = Dict{Tuple{T, T}, EdgeType.EType}()
    for e in edges(nzxd.mg)
        et[(src(e), dst(e))] = EdgeType.SIM
    end
    zxg = ZXGraph{T, P}(nzxd.mg, nzxd.ps, nzxd.st, et, nzxd.layout, nzxd.phase_ids, nzxd.scalar, zxd, nzxd.inputs, nzxd.outputs)

    for e in eH
        v1, v2 = e
        add_edge!(zxg, v1, v2)
    end

    return zxg
end

LightGraphs.has_edge(zxg::ZXGraph, vs...) = has_edge(zxg.mg, vs...)
LightGraphs.nv(zxg::ZXGraph) = nv(zxg.mg)
LightGraphs.ne(zxg::ZXGraph) = ne(zxg.mg)
LightGraphs.outneighbors(zxg::ZXGraph, v::Integer) = outneighbors(zxg.mg, v)
LightGraphs.inneighbors(zxg::ZXGraph, v::Integer) = inneighbors(zxg.mg, v)
LightGraphs.neighbors(zxg::ZXGraph, v::Integer) = neighbors(zxg.mg, v)
LightGraphs.degree(zxg::ZXGraph, v::Integer) = degree(zxg.mg, v)
LightGraphs.indegree(zxg::ZXGraph, v::Integer) = degree(zxg, v)
LightGraphs.outdegree(zxg::ZXGraph, v::Integer) = degree(zxg, v)
function LightGraphs.rem_edge!(zxg::ZXGraph, v1::Integer, v2::Integer)
    if rem_edge!(zxg.mg, v1, v2)
        delete!(zxg.et, (min(v1, v2), max(v1, v2)))
        return true
    end
    return false
end

function LightGraphs.add_edge!(zxg::ZXGraph, v1::Integer, v2::Integer, edge_type::EdgeType.EType = EdgeType.HAD)
    if has_vertex(zxg.mg, v1) && has_vertex(zxg.mg, v2)
        if v1 == v2
            if edge_type == EdgeType.HAD
                set_phase!(zxg, v1, phase(zxg, v1)+1)
                add_power!(zxg, -1)
            end
            return true
        else
            if has_edge(zxg, v1, v2)
                if is_hadamard(zxg, v1, v2)
                    add_power!(zxg, -2)
                    return rem_edge!(zxg, v1, v2)
                else
                    return false
                end
            elseif add_edge!(zxg.mg, v1, v2)
                zxg.et[(min(v1, v2), max(v1, v2))] = edge_type
                return true
            end
        end
    end
    return false
end

spider_type(zxg::ZXGraph, v::Integer) = zxg.st[v]
phase(zxg::ZXGraph, v::Integer) = zxg.ps[v]
function set_phase!(zxg::ZXGraph{T, P}, v::T, p::P) where {T, P}
    if has_vertex(zxg.mg, v)
        while p < 0
            p += 2
        end
        p = rem(p, 2)
        zxg.ps[v] = p
        return true
    end
    return false
end
nqubits(zxg::ZXGraph) = zxg.layout.nbits

qubit_loc(zxg::ZXGraph{T, P}, v::T) where {T, P} = qubit_loc(zxg.layout, v)
function column_loc(zxg::ZXGraph{T, P}, v::T) where {T, P}
    c_loc = column_loc(zxg.layout, v)
    if c_loc !== nothing
        if spider_type(zxg, v) == SpiderType.Out
            nb = neighbors(zxg, v)
            if length(nb) == 1
                nb = nb[1]
                spider_type(zxg, nb) == SpiderType.In && return 3//1
                c_loc = floor(column_loc(zxg, nb) + 2)
            else
                c_loc = 1000
            end
        end
        if spider_type(zxg, v) == SpiderType.In
            nb = neighbors(zxg, v)[1]
            spider_type(zxg, nb) == SpiderType.Out && return 1//1
            c_loc = ceil(column_loc(zxg, nb) - 2)
        end
    end
    return c_loc
end

function is_hadamard(zxg::ZXGraph, v1::Integer, v2::Integer)
    if has_edge(zxg, v1, v2)
        src = min(v1, v2)
        dst = max(v1, v2)
        return zxg.et[(src, dst)] == EdgeType.HAD
    end
    return false
end
spiders(zxg::ZXGraph) = vertices(zxg.mg)

function rem_spiders!(zxg::ZXGraph{T, P}, vs::Vector{T}) where {T, P}
    if rem_vertices!(zxg.mg, vs)
        for v in vs
            delete!(zxg.ps, v)
            delete!(zxg.st, v)
            delete!(zxg.phase_ids, v)
            rem_vertex!(zxg.layout, v)
        end
        return true
    end
    return false
end
rem_spider!(zxg::ZXGraph{T, P}, v::T) where {T, P} = rem_spiders!(zxg, [v])

function add_spider!(zxg::ZXGraph{T, P}, st::SpiderType.SType, phase::P = zero(P), connect::Vector{T}=T[]) where {T<:Integer, P}
    v = add_vertex!(zxg.mg)[1]
    set_phase!(zxg, v, phase)
    zxg.st[v] = st
    if st in (SpiderType.Z, SpiderType.X)
        zxg.phase_ids[v] = (v, 1)
    end
    if all(has_vertex(zxg.mg, c) for c in connect)
        for c in connect
            add_edge!(zxg, v, c)
        end
    end
    return v
end
function insert_spider!(zxg::ZXGraph{T, P}, v1::T, v2::T, phase::P = zero(P)) where {T<:Integer, P}
    v = add_spider!(zxg, SpiderType.Z, phase, [v1, v2])
    rem_edge!(zxg, v1, v2)
    return v
end

function print_spider(io::IO, zxg::ZXGraph{T}, v::T) where {T<:Integer}
    st_v = spider_type(zxg, v)
    if st_v == SpiderType.Z
        printstyled(io, "S_$(v){phase = $(phase(zxg, v))"*(zxg.ps[v] isa Phase ? "}" : "⋅π}"); color = :green)
    elseif st_v == SpiderType.In
        print(io, "S_$(v){input}")
    elseif st_v == SpiderType.Out
        print(io, "S_$(v){output}")
    end
end

function Base.show(io::IO, zxg::ZXGraph{T}) where {T<:Integer}
    println(io, "ZX-graph with $(nv(zxg)) vertices and $(ne(zxg)) edges:")
    vs = sort!(spiders(zxg))
    for i = 1:length(vs)
        for j = i+1:length(vs)
            if has_edge(zxg, vs[i], vs[j])
                print(io, "(")
                print_spider(io, zxg, vs[i])
                if is_hadamard(zxg, vs[i], vs[j])
                    printstyled(io, " <-> "; color = :blue)
                else
                    print(io, " <-> ")
                end
                print_spider(io, zxg, vs[j])
                print(io, ")\n")
            end
        end
    end
end

function round_phases!(zxg::ZXGraph{T, P}) where {T<:Integer, P}
    ps = zxg.ps
    for v in keys(ps)
        while ps[v] < 0
            ps[v] += 2
        end
        ps[v] = rem(ps[v], 2)
    end
end

"""
    is_interior(zxg::ZXGraph, v)

Return `true` if `v` is a interior spider of `zxg`.
"""
function is_interior(zxg::ZXGraph{T, P}, v::T) where {T, P}
    if has_vertex(zxg.mg, v)
        (spider_type(zxg, v) == SpiderType.In || spider_type(zxg, v) == SpiderType.Out) && return false
        for u in neighbors(zxg, v)
            if spider_type(zxg, u) == SpiderType.In || spider_type(zxg, u) == SpiderType.Out
                return false
            end
        end
        return true
    end
    return false
end

get_outputs(zxg::ZXGraph) = get_outputs(zxg.master)
get_inputs(zxg::ZXGraph) = get_inputs(zxg.master)

# TODO: remove it?
function spider_sequence(zxg::ZXGraph{T, P}) where {T, P}
    nbits = nqubits(zxg)
    if nbits > 0
        vs = spiders(zxg)
        spider_seq = Vector{Vector{T}}(undef, nbits)
        for q = 1:nbits
            spider_seq[q] = Vector{T}()
        end
        for v in vs
            if qubit_loc(zxg, v) !== nothing
                q_loc = Int(qubit_loc(zxg, v))
                q_loc > 0 && push!(spider_seq[q_loc], v)
            end
        end
        for q = 1:nbits
            sort!(spider_seq[q], by = (v -> column_loc(zxg, v)))
        end
        return spider_seq
    end
end

function generate_layout!(zxg::ZXGraph{T, P}) where {T, P}
    layout = zxg.layout
    nbits = length(zxg.inputs)
    vs_frontier = copy(zxg.inputs)
    vs_generated = Set(vs_frontier)
    frontier_col = [1//1 for _ = 1:nbits]
    frontier_active = [true for _ = 1:nbits]
    for i = 1:nbits
        set_qubit!(layout, vs_frontier[i], i)
        set_column!(layout, vs_frontier[i], 1//1)
    end
    gad_col = 1//1
    for v in spiders(zxg)
        if degree(zxg, v) == 1 && spider_type(zxg, v) == SpiderType.Z
            v1 = neighbors(zxg, v)[1]
            set_loc!(layout, v, -1//1, gad_col)
            set_loc!(layout, v1, 0//1, gad_col)
            push!(vs_generated, v, v1)
            gad_col += 1
        end
    end

    while !(zxg.outputs ⊆ vs_frontier)
        while any(frontier_active)
            for q in 1:nbits
                if frontier_active[q]
                    v = vs_frontier[q]
                    nb = neighbors(zxg, v)
                    nb_not_gen = nb[findall(v -> !(v in vs_generated), nb)]
                    if length(nb_not_gen) == 1
                        set_loc!(layout, v, q, frontier_col[q])
                        frontier_col[q] += 1
                        push!(vs_generated, v)
                        vs_frontier[q] = nb_not_gen[1]
                    else
                        frontier_active[q] = false
                    end
                end
            end
        end
        for q = 1:nbits
            isupdated = false
            v = vs_frontier[q]
            nb = neighbors(zxg, v)
            for v1 in nb
                if !(v1 in vs_generated)
                    q1 = findfirst(isequal(v1), vs_frontier)
                    if q1 !== nothing
                        col = maximum(frontier_col[min(q, q1):max(q, q1)])
                        set_loc!(layout, v, q, col)
                        set_loc!(layout, v1, q1, col)
                        nb_v1 = neighbors(zxg, v1)
                        push!(vs_generated, v, v1)
                        new_v1 = nb_v1[findall(v -> !(v in vs_generated), nb_v1)]
                        new_v = nb[findall(v -> !(v in vs_generated), nb)]
                        if length(new_v) == 1
                            vs_frontier[q] = new_v[]
                            frontier_active[q] = true
                            isupdated = true
                            for i in min(q, q1):max(q, q1)
                                frontier_col[i] = col + 1
                            end
                        else
                            delete!(vs_generated, v)
                        end
                        if length(new_v1) == 1
                            vs_frontier[q1] = new_v1[]
                            frontier_active[q1] = true
                            isupdated = true
                            for i in min(q, q1):max(q, q1)
                                frontier_col[i] = col + 1
                            end
                        else
                            delete!(vs_generated, v1)
                        end
                        break
                    end
                end
            end
            isupdated && break
        end
    end
    for q = 1:length(vs_frontier)
        set_loc!(layout, vs_frontier[q], q, frontier_col[q])
    end
    return layout
end

scalar(zxg::ZXGraph) = zxg.scalar

function add_global_phase!(zxg::ZXGraph{T, P}, p::P) where {T, P}
    add_phase!(zxg.scalar, p)
    return zxg
end

function add_power!(zxg::ZXGraph, n)
    add_power!(zxg.scalar, n)
    return zxg
end
