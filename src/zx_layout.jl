"""
    ZXLayout

A struct for the layout information of `ZXDiagram` and `ZXGraph`.
"""
struct ZXLayout{T<:Integer}
    nbits::Int
    spider_q::Dict{T, Rational{Int}}
    spider_col::Dict{T, Rational{Int}}
end

ZXLayout(nbits::Integer, spider_q::Dict{T, Rational{Int}}, spider_col::Dict{T, Rational{Int}}) where {T} = ZXLayout{T}(Int(nbits), spider_q, spider_col)
ZXLayout{T}() where {T} = ZXLayout(0, Dict{T, Rational{Int}}(), Dict{T, Rational{Int}}())

Base.copy(layout::ZXLayout) = ZXLayout(layout.nbits, copy(layout.spider_q), copy(layout.spider_col))
function LightGraphs.rem_vertex!(layout::ZXLayout{T}, v::T) where {T}
    delete!(layout.spider_q, v)
    delete!(layout.spider_col, v)
    return
end

"""
    qubit_loc(layout, v)

Return the qubit number corresponding to the spider `v`.
"""
qubit_loc(layout::ZXLayout{T}, v::T) where T = get(layout.spider_q, v, nothing)

"""
    column_loc(layout, v)

Return the column number corresponding to the spider `v`.
"""
column_loc(layout::ZXLayout{T}, v::T) where T = get(layout.spider_col, v, nothing)

"""
    set_qubit!(layout, v, q)

Set the qubit number of the spider `v`.
"""
function set_qubit!(layout::ZXLayout{T}, v::T, q) where T
    layout.spider_q[v] = q
    return layout
end

"""
    set_qubit!(layout, v, q)

Set the column number of the spider `v`.
"""
function set_column!(layout::ZXLayout{T}, v::T, col) where T
    layout.spider_col[v] = col
    return layout
end

"""
    set_loc!(layout, v, q, col)

Set the location of the spider `v`.
"""
function set_loc!(layout::ZXLayout{T}, v::T, q, col) where T
    set_qubit!(layout, v, q)
    set_column!(layout, v, col)
    return layout
end
