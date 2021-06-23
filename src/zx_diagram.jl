module SpiderType
    @enum SType Z X H In Out
end  # module SpiderType

"""
    ZXDiagram{T, P}
This is the type for representing ZX-diagrams.
"""
struct ZXDiagram{T<:Integer, P} <: AbstractZXDiagram{T, P}
    mg::Multigraph{T}

    st::Dict{T, SpiderType.SType}
    ps::Dict{T, P}

    layout::ZXLayout{T}
    phase_ids::Dict{T, Tuple{T, Int}}

    scalar::Scalar{P}
    inputs::Vector{T}
    outputs::Vector{T}

    function ZXDiagram{T, P}(mg::Multigraph{T}, st::Dict{T, SpiderType.SType}, ps::Dict{T, P},
        layout::ZXLayout{T}, phase_ids::Dict{T, Tuple{T, Int}} = Dict{T, Tuple{T, Int}}(),
        s::Scalar{P} = Scalar{P}(), 
        inputs::Vector{T} = Vector{T}(), outputs::Vector{T} = Vector{T}()) where {T<:Integer, P}
        if nv(mg) == length(ps) && nv(mg) == length(st)
            if length(phase_ids) == 0
                for v in vertices(mg)
                    if st[v] in (SpiderType.Z, SpiderType.X)
                        phase_ids[v] = (v, 1)
                    end
                end
            end
            if length(inputs) == 0
                for v in vertices(mg)
                    if st[v] == SpiderType.In
                        push!(inputs, v)
                    end
                end
                if layout.nbits > 0
                    sort!(inputs, by = (v -> qubit_loc(layout, v)))
                end
            end
            if length(outputs) == 0
                for v in vertices(mg)
                    if st[v] == SpiderType.Out
                        push!(outputs, v)
                    end
                end
                if layout.nbits > 0
                    sort!(outputs, by = (v -> qubit_loc(layout, v)))
                end
            end
            zxd = new{T, P}(mg, st, ps, layout, phase_ids, s, inputs, outputs)
            round_phases!(zxd)
            return zxd
        else
            error("There should be a phase and a type for each spider!")
        end
    end
end

"""
    ZXDiagram(mg::Multigraph{T}, st::Dict{T, SpiderType.SType}, ps::Dict{T, P},
        layout::ZXLayout{T} = ZXLayout{T}(),
        phase_ids::Dict{T,Tuple{T, Int}} = Dict{T,Tuple{T,Int}}()) where {T, P}
    ZXDiagram(mg::Multigraph{T}, st::Vector{SpiderType.SType}, ps::Vector{P},
        layout::ZXLayout{T} = ZXLayout{T}()) where {T, P}

Construct a ZXDiagram with all information.

```jldoctest
julia> using LightGraphs, Multigraphs, ZXCalculus;

julia> using ZXCalculus.SpiderType: In, Out, H, Z, X;

julia> mg = Multigraph(5);

julia> for i = 1:4
           add_edge!(mg, i, i+1)
       end;

julia> ZXDiagram(mg, [In, Z, H, X, Out], [0//1, 1, 0, 1//2, 0])
ZX-diagram with 5 vertices and 4 multiple edges:
(S_1{input} <-1-> S_2{phase = 1//1⋅π})
(S_2{phase = 1//1⋅π} <-1-> S_3{H})
(S_3{H} <-1-> S_4{phase = 1//2⋅π})
(S_4{phase = 1//2⋅π} <-1-> S_5{output})

```
"""
ZXDiagram(mg::Multigraph{T}, st::Dict{T, SpiderType.SType}, ps::Dict{T, P},
    layout::ZXLayout{T} = ZXLayout{T}(),
    phase_ids::Dict{T,Tuple{T, Int}} = Dict{T,Tuple{T,Int}}()) where {T, P} = ZXDiagram{T, P}(mg, st, ps, layout, phase_ids)
ZXDiagram(mg::Multigraph{T}, st::Vector{SpiderType.SType}, ps::Vector{P},
    layout::ZXLayout{T} = ZXLayout{T}()) where {T, P} =
    ZXDiagram(mg, Dict(zip(sort!(vertices(mg)), st)), Dict(zip(sort!(vertices(mg)), ps)), layout)

"""
    ZXDiagram(nbits)

Construct a ZXDiagram of a empty circuit with qubit number `nbit`

```jldoctest; setup = :(using ZXCalculus)
julia> zxd = ZXDiagram(3)
ZX-diagram with 6 vertices and 3 multiple edges:
(S_1{input} <-1-> S_2{output})
(S_3{input} <-1-> S_4{output})
(S_5{input} <-1-> S_6{output})

```
"""
function ZXDiagram(nbits::T) where {T<:Integer}
    mg = Multigraph(2*nbits)
    st = [SpiderType.In for _ = 1:2*nbits]
    ps = [Phase(0//1) for _ = 1:2*nbits]
    spider_q = Dict{T, Rational{Int}}()
    spider_col = Dict{T, Rational{Int}}()
    for i = 1:nbits
        add_edge!(mg, 2*i-1, 2*i)
        @inbounds st[2*i] = SpiderType.Out
        spider_q[2*i-1] = i
        spider_col[2*i-1] = 2
        spider_q[2*i] = i
        spider_col[2*i] = 1
    end
    layout = ZXLayout(nbits, spider_q, spider_col)
    return ZXDiagram(mg, st, ps, layout)
end

Base.copy(zxd::ZXDiagram{T, P}) where {T, P} = ZXDiagram{T, P}(copy(zxd.mg), copy(zxd.st), copy(zxd.ps), copy(zxd.layout),
    deepcopy(zxd.phase_ids), copy(zxd.scalar), copy(zxd.inputs), copy(zxd.outputs))

"""
    spider_type(zxd, v)

Returns the spider type of a spider.
"""
spider_type(zxd::ZXDiagram{T, P}, v::T) where {T<:Integer, P} = zxd.st[v]

"""
    phase(zxd, v)

Returns the phase of a spider. If the spider is not a Z or X spider, then return 0.
"""
phase(zxd::ZXDiagram{T, P}, v::T) where {T<:Integer, P} = zxd.ps[v]


"""
    set_phase!(zxd, v, p)

Set the phase of `v` in `zxd` to `p`.
"""
function set_phase!(zxd::ZXDiagram{T, P}, v::T, p::P) where {T, P}
    if has_vertex(zxd.mg, v)
        while p < 0
            p += 2
        end
        p = rem(p, 2)
        zxd.ps[v] = p
        return true
    end
    return false
end


"""
    nqubits(zxd)

Returns the qubit number of a ZX-diagram.
"""
nqubits(zxd::ZXDiagram) = zxd.layout.nbits

function print_spider(io::IO, zxd::ZXDiagram{T, P}, v::T) where {T<:Integer, P}
    st_v = spider_type(zxd, v)
    if st_v == SpiderType.Z
        printstyled(io, "S_$(v){phase = $(zxd.ps[v])"*(zxd.ps[v] isa Phase ? "}" : "⋅π}"); color = :green)
    elseif st_v == SpiderType.X
        printstyled(io, "S_$(v){phase = $(zxd.ps[v])"*(zxd.ps[v] isa Phase ? "}" : "⋅π}"); color = :red)
    elseif st_v == SpiderType.H
        printstyled(io, "S_$(v){H}"; color = :yellow)
    elseif st_v == SpiderType.In
        print(io, "S_$(v){input}")
    elseif st_v == SpiderType.Out
        print(io, "S_$(v){output}")
    end
end

function Base.show(io::IO, zxd::ZXDiagram{T, P}) where {T<:Integer, P}
    println(io, "ZX-diagram with $(nv(zxd.mg)) vertices and $(ne(zxd.mg)) multiple edges:")
    for v1 in sort!(vertices(zxd.mg))
        for v2 in neighbors(zxd.mg, v1)
            if v2 >= v1
                print(io, "(")
                print_spider(io, zxd, v1)
                print(io, " <-$(mul(zxd.mg, v1, v2))-> ")
                print_spider(io, zxd, v2)
                print(io, ")\n")
            end
        end
    end
end

"""
    nv(zxd)

Returns the number of vertices (spiders) of a ZX-diagram.
"""
LightGraphs.nv(zxd::ZXDiagram) = nv(zxd.mg)

"""
    ne(zxd; count_mul = false)

Returns the number of edges of a ZX-diagram. If `count_mul`, it will return the
sum of multiplicities of all multiple edges. Otherwise, it will return the
number of multiple edges.
"""
LightGraphs.ne(zxd::ZXDiagram; count_mul::Bool = false) = ne(zxd.mg, count_mul = count_mul)

LightGraphs.outneighbors(zxd::ZXDiagram, v; count_mul::Bool = false) = outneighbors(zxd.mg, v, count_mul = count_mul)
LightGraphs.inneighbors(zxd::ZXDiagram, v; count_mul::Bool = false) = inneighbors(zxd.mg, v, count_mul = count_mul)

LightGraphs.degree(zxd::ZXDiagram, v::Integer) = degree(zxd.mg, v)
LightGraphs.indegree(zxd::ZXDiagram, v::Integer) = degree(zxd, v)
LightGraphs.outdegree(zxd::ZXDiagram, v::Integer) = degree(zxd, v)

"""
    neighbors(zxd, v; count_mul = false)

Returns a vector of vertices connected to `v`. If `count_mul`, there will be
multiple copy for each vertex. Otherwise, each vertex will only appear once.
"""
LightGraphs.neighbors(zxd::ZXDiagram, v; count_mul::Bool = false) = neighbors(zxd.mg, v, count_mul = count_mul)
function LightGraphs.rem_edge!(zxd::ZXDiagram, x...)
    rem_edge!(zxd.mg, x...)
end
function LightGraphs.add_edge!(zxd::ZXDiagram, x...)
    add_edge!(zxd.mg, x...)
end

"""
    rem_spiders!(zxd, vs)

Remove spiders indexed by `vs`.
"""
function rem_spiders!(zxd::ZXDiagram{T, P}, vs::Vector{T}) where {T<:Integer, P}
    if rem_vertices!(zxd.mg, vs)
        for v in vs
            delete!(zxd.ps, v)
            delete!(zxd.st, v)
            delete!(zxd.phase_ids, v)
            rem_vertex!(zxd.layout, v)
        end
        return true
    end
    return false
end

"""
    rem_spider!(zxd, v)

Remove a spider indexed by `v`.
"""
rem_spider!(zxd::ZXDiagram{T, P}, v::T) where {T<:Integer, P} = rem_spiders!(zxd, [v])

"""
    add_spider!(zxd, spider_type, phase = 0, connect = [])

Add a new spider which is of the type `spider_type` with phase `phase` and
connected to the vertices `connect`.
"""
function add_spider!(zxd::ZXDiagram{T, P}, st::SpiderType.SType, phase::P = zero(P), connect::Vector{T}=T[]) where {T<:Integer, P}
    v = add_vertex!(zxd.mg)[1]
    set_phase!(zxd, v, phase)
    zxd.st[v] = st
    if st in (SpiderType.Z, SpiderType.X)
        zxd.phase_ids[v] = (v, 1)
    end
    if all(has_vertex(zxd.mg, c) for c in connect)
        for c in connect
            add_edge!(zxd.mg, v, c)
        end
    end
    return v
end

"""
    insert_spider!(zxd, v1, v2, spider_type, phase = 0)

Insert a spider of the type `spider_type` with phase = `phase`, between two
vertices `v1` and `v2`. It will insert multiple times if the edge between
`v1` and `v2` is a multiple edge. Also it will remove the original edge between
`v1` and `v2`.
"""
function insert_spider!(zxd::ZXDiagram{T, P}, v1::T, v2::T, st::SpiderType.SType, phase::P = zero(P)) where {T<:Integer, P}
    mt = mul(zxd.mg, v1, v2)
    vs = Vector{T}(undef, mt)
    for i = 1:mt
        v = add_spider!(zxd, st, phase, [v1, v2])
        @inbounds vs[i] = v
        rem_edge!(zxd, v1, v2)
    end
    return vs
end

"""
    round_phases!(zxd)

Round phases between [0, 2π).
"""
function round_phases!(zxd::ZXDiagram{T, P}) where {T<:Integer, P}
    ps = zxd.ps
    for v in keys(ps)
        while ps[v] < 0
            ps[v] += 2
        end
        ps[v] = rem(ps[v], 2)
    end
    return
end

spiders(zxd::ZXDiagram) = vertices(zxd.mg)
qubit_loc(zxd::ZXDiagram{T, P}, v::T) where {T, P} = qubit_loc(zxd.layout, v)
function column_loc(zxd::ZXDiagram{T, P}, v::T) where {T, P}
    c_loc = column_loc(zxd.layout, v)
    return c_loc
end

"""
    push_gate!(zxd, ::Val{M}, locs...[, phase]; autoconvert=true)

Push an `M` gate to the end of qubit `loc` where `M` can be `:Z`, `:X`, `:H`, `:SWAP`, `:CNOT` and `:CZ`.
If `M` is `:Z` or `:X`, `phase` will be available and it will push a
rotation `M` gate with angle `phase * π`.
If `autoconvert` is `false`, the input `phase` should be a rational numbers.
"""
function push_gate!(zxd::ZXDiagram{T, P}, ::Val{:Z}, loc::T, phase = zero(P); autoconvert::Bool=true) where {T, P}
    @inbounds out_id = get_outputs(zxd)[loc]
    @inbounds bound_id = neighbors(zxd, out_id)[1]
    rphase = autoconvert ? safe_convert(P, phase) : phase
    insert_spider!(zxd, bound_id, out_id, SpiderType.Z, rphase)
    return zxd
end

function push_gate!(zxd::ZXDiagram{T, P}, ::Val{:X}, loc::T, phase = zero(P); autoconvert::Bool=true) where {T, P}
    @inbounds out_id = get_outputs(zxd)[loc]
    @inbounds bound_id = neighbors(zxd, out_id)[1]
    rphase = autoconvert ? safe_convert(P, phase) : phase
    insert_spider!(zxd, bound_id, out_id, SpiderType.X, rphase)
    return zxd
end

function push_gate!(zxd::ZXDiagram{T, P}, ::Val{:H}, loc::T) where {T, P}
    @inbounds out_id = get_outputs(zxd)[loc]
    @inbounds bound_id = neighbors(zxd, out_id)[1]
    insert_spider!(zxd, bound_id, out_id, SpiderType.H)
    return zxd
end

function push_gate!(zxd::ZXDiagram{T, P}, ::Val{:SWAP}, locs::Vector{T}) where {T, P}
    q1, q2 = locs
    push_gate!(zxd, Val{:Z}(), q1)
    push_gate!(zxd, Val{:Z}(), q2)
    push_gate!(zxd, Val{:Z}(), q1)
    push_gate!(zxd, Val{:Z}(), q2)
    v1, v2, bound_id1, bound_id2 = (sort!(spiders(zxd)))[end-3:end]
    rem_edge!(zxd, v1, bound_id1)
    rem_edge!(zxd, v2, bound_id2)
    add_edge!(zxd, v1, bound_id2)
    add_edge!(zxd, v2, bound_id1)
    return zxd
end

function push_gate!(zxd::ZXDiagram{T, P}, ::Val{:CNOT}, loc::T, ctrl::T) where {T, P}
    push_gate!(zxd, Val{:Z}(), ctrl)
    push_gate!(zxd, Val{:X}(), loc)
    @inbounds v1, v2 = (sort!(spiders(zxd)))[end-1:end]
    add_edge!(zxd, v1, v2)
    add_power!(zxd, 1)
    return zxd
end

function push_gate!(zxd::ZXDiagram{T, P}, ::Val{:CZ}, loc::T, ctrl::T) where {T, P}
    push_gate!(zxd, Val{:Z}(), ctrl)
    push_gate!(zxd, Val{:Z}(), loc)
    @inbounds v1, v2 = (sort!(spiders(zxd)))[end-1:end]
    add_edge!(zxd, v1, v2)
    insert_spider!(zxd, v1, v2, SpiderType.H)
    add_power!(zxd, 1)
    return zxd
end

"""
    pushfirst_gate!(zxd, ::Val{M}, loc[, phase])

Push an `M` gate to the beginning of qubit `loc` where `M` can be `:Z`, `:X`, `:H`, `:SWAP`, `:CNOT` and `:CZ`.
If `M` is `:Z` or `:X`, `phase` will be available and it will push a
rotation `M` gate with angle `phase * π`.
"""
function pushfirst_gate!(zxd::ZXDiagram{T, P}, ::Val{:Z}, loc::T, phase::P = zero(P)) where {T, P}
    @inbounds in_id = get_inputs(zxd)[loc]
    @inbounds bound_id = neighbors(zxd, in_id)[1]
    insert_spider!(zxd, in_id, bound_id, SpiderType.Z, phase)
    return zxd
end

function pushfirst_gate!(zxd::ZXDiagram{T, P}, ::Val{:X}, loc::T, phase::P = zero(P)) where {T, P}
    @inbounds in_id = get_inputs(zxd)[loc]
    @inbounds bound_id = neighbors(zxd, in_id)[1]
    insert_spider!(zxd, in_id, bound_id, SpiderType.X, phase)
    return zxd
end

function pushfirst_gate!(zxd::ZXDiagram{T, P}, ::Val{:H}, loc::T) where {T, P}
    @inbounds in_id = get_inputs(zxd)[loc]
    @inbounds bound_id = neighbors(zxd, in_id)[1]
    insert_spider!(zxd, in_id, bound_id, SpiderType.H)
    return zxd
end

function pushfirst_gate!(zxd::ZXDiagram{T, P}, ::Val{:SWAP}, locs::Vector{T}) where {T, P}
    q1, q2 = locs
    pushfirst_gate!(zxd, Val{:Z}(), q1)
    pushfirst_gate!(zxd, Val{:Z}(), q2)
    pushfirst_gate!(zxd, Val{:Z}(), q1)
    pushfirst_gate!(zxd, Val{:Z}(), q2)
    @inbounds v1, v2, bound_id1, bound_id2 = (sort!(spiders(zxd)))[end-3:end]
    rem_edge!(zxd, v1, bound_id1)
    rem_edge!(zxd, v2, bound_id2)
    add_edge!(zxd, v1, bound_id2)
    add_edge!(zxd, v2, bound_id1)
    return zxd
end

function pushfirst_gate!(zxd::ZXDiagram{T, P}, ::Val{:CNOT}, loc::T, ctrl::T) where {T, P}
    pushfirst_gate!(zxd, Val{:Z}(), ctrl)
    pushfirst_gate!(zxd, Val{:X}(), loc)
    @inbounds v1, v2 = (sort!(spiders(zxd)))[end-1:end]
    add_edge!(zxd, v1, v2)
    add_power!(zxd, 1)
    return zxd
end

function pushfirst_gate!(zxd::ZXDiagram{T, P}, ::Val{:CZ}, loc::T, ctrl::T) where {T, P}
    pushfirst_gate!(zxd, Val{:Z}(), ctrl)
    pushfirst_gate!(zxd, Val{:Z}(), loc)
    @inbounds v1, v2 = (sort!(spiders(zxd)))[end-1:end]
    add_edge!(zxd, v1, v2)
    insert_spider!(zxd, v1, v2, SpiderType.H)
    add_power!(zxd, 1)
    return zxd
end

"""
    tcount(zxd)

Returns the T-count of a ZX-diagram.
"""
tcount(cir::AbstractZXDiagram) = sum([phase(cir, v) % 1//2 != 0 for v in spiders(cir)])

"""
    get_inputs(zxd)

Returns a vector of input ids.
"""
get_inputs(zxd::ZXDiagram) = zxd.inputs

"""
    get_outputs(zxd)

Returns a vector of output ids.
"""
get_outputs(zxd::ZXDiagram) = zxd.outputs

"""
    scalar(zxd)

Returns the scalar of `zxd`.
"""
scalar(zxd::ZXDiagram) = zxd.scalar

function add_global_phase!(zxd::ZXDiagram{T, P}, p::P) where {T, P}
    add_phase!(zxd.scalar, p)
    return zxd
end

function add_power!(zxd::ZXDiagram, n)
    add_power!(zxd.scalar, n)
    return zxd
end

function spider_sequence(zxd::ZXDiagram{T, P}) where {T, P}
    seq = []
    generate_layout!(zxd, seq)
    return seq
end

function generate_layout!(zxd::ZXDiagram{T, P}, seq::Vector{Any} = []) where {T, P}
    layout = zxd.layout
    nbits = length(zxd.inputs)
    vs_frontier = copy(zxd.inputs)
    vs_generated = Set(vs_frontier)
    frontier_col = [1//1 for _ = 1:nbits]
    frontier_active = [true for _ = 1:nbits]
    for i = 1:nbits
        set_qubit!(layout, vs_frontier[i], i)
        set_column!(layout, vs_frontier[i], 1//1)
    end

    while !(zxd.outputs ⊆ vs_frontier)
        while any(frontier_active)
            for q in 1:nbits
                if frontier_active[q]
                    v = vs_frontier[q]
                    nb = neighbors(zxd, v)
                    if length(nb) <= 2
                        set_loc!(layout, v, q, frontier_col[q])
                        push!(seq, v)
                        push!(vs_generated, v)
                        q_active = false
                        for v1 in nb
                            if !(v1 in vs_generated)
                                vs_frontier[q] = v1
                                frontier_col[q] += 1
                                q_active = true
                                break
                            end
                        end
                        frontier_active[q] = q_active
                    else
                        frontier_active[q] = false
                    end
                end
            end
        end
        for q = 1:nbits
            v = vs_frontier[q]
            nb = neighbors(zxd, v)
            isupdated = false
            for v1 in nb
                if !(v1 in vs_generated)
                    q1 = findfirst(isequal(v1), vs_frontier)
                    if q1 !== nothing
                        col = maximum(frontier_col[min(q, q1):max(q, q1)])
                        set_loc!(layout, v, q, col)
                        set_loc!(layout, v1, q1, col)
                        push!(vs_generated, v, v1)
                        push!(seq, (v, v1))
                        nb_v1 = neighbors(zxd, v1)
                        new_v1 = nb_v1[findfirst(v -> !(v in vs_generated), nb_v1)]
                        new_v = nb[findfirst(v -> !(v in vs_generated), nb)]
                        vs_frontier[q] = new_v
                        vs_frontier[q1] = new_v1
                        for i in min(q, q1):max(q, q1)
                            frontier_col[i] = col + 1
                        end
                        frontier_active[q] = true
                        frontier_active[q1] = true
                        isupdated = true
                        break
                    elseif spider_type(zxd, v1) == SpiderType.H && degree(zxd, v1) == 2
                        nb_v1 = neighbors(zxd, v1)
                        v2 = nb_v1[findfirst(!isequal(v), nb_v1)]
                        q2 = findfirst(isequal(v2), vs_frontier)
                        if q2 !== nothing
                            col = maximum(frontier_col[min(q, q2):max(q, q2)])
                            set_loc!(layout, v, q, col)
                            set_loc!(layout, v2, q2, col)
                            q1 = (q + q2)//2
                            denominator(q1) == 1 && (q1 += 1//2)
                            set_loc!(layout, v1, q1, col)
                            push!(vs_generated, v, v1, v2)
                            push!(seq, (v, v1, v2))
                            nb_v2 = neighbors(zxd, v2)
                            new_v = nb[findfirst(v -> !(v in vs_generated), nb)]
                            new_v2 = nb_v2[findfirst(v -> !(v in vs_generated), nb_v2)]
                            vs_frontier[q] = new_v
                            vs_frontier[q2] = new_v2
                            for i in min(q, q2):max(q, q2)
                                frontier_col[i] = col + 1
                            end
                            frontier_active[q] = true
                            frontier_active[q2] = true
                            isupdated = true
                            break
                        end
                    end
                end
                isupdated && break
            end
        end
    end
    for q = 1:length(vs_frontier)
        set_loc!(layout, vs_frontier[q], q, frontier_col[q])
    end
    return layout
end

"""
    continued_fraction(ϕ, n::Int) -> Rational

Obtain `s` and `r` from `ϕ` that satisfies `|s//r - ϕ| ≦ 1/2r²`
"""
function continued_fraction(fl, n::Int)
    if n == 1 || abs(mod(fl, 1)) < 1e-10
        Rational(floor(Int, fl), 1)
    else
        floor(Int, fl) + 1//continued_fraction(1/mod(fl, 1), n-1)
    end
end

safe_convert(::Type{T}, x) where T = convert(T, x)
safe_convert(::Type{T}, x::T) where T<:Rational = x
function safe_convert(::Type{T}, x::Real) where T<:Rational
    local fr
    for n=1:16 # at most 20 steps, otherwise the number may overflow.
        fr = continued_fraction(x, n)
        abs(fr - x) < 1e-12 && return fr
    end
    @warn "converting phase to rational, but with rounding error $(abs(fr-x))."
    return fr
end
