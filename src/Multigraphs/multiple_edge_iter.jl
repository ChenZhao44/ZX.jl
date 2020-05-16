using LightGraphs

import Base: eltype, iterate, length

export MultipleEdgeIter

struct MultipleEdgeIter{G<:AbstractMultigraph} <:AbstractEdgeIter
    g::G
end

eltype(::Type{MultipleEdgeIter{AbstractMultigraph{T}}}) where {T<:Integer} = MultipleEdge{T, Int}

function iterate(eit::MultipleEdgeIter{G}, state=(one(eltype(eit.g)), one(eltype(eit.g)))) where {G <: AbstractMultigraph}
    g = eit.g
    n = nv(g)
    u, i = state

    @inbounds while u <= n
        list_u = outneighbors(g, u)
        if i > length(list_u)
            u += 1
            i = one(u)
            continue
        end
        if is_directed(g)
            e = MultipleEdge(u, list_u[i], g.adjmx[u, list_u[i]])
            state = (u, i + 1)
            return e, state
        else
            if list_u[i] >= u
                e = MultipleEdge(u, list_u[i], g.adjmx[list_u[i], u])
                state = (u, i + 1)
                return e, state
            else
                i += 1
            end
        end
    end

    if n == 0 || u > n
        return nothing
    end
end

function iterate(eit::MultipleEdgeIter{G}, state=(one(eltype(eit.g)), one(eltype(eit.g)))) where {G <: Multigraph}
    g = eit.g
    n = nv(g)
    vs = vertices(g)
    u, i = state

    @inbounds while u <= n
        list_u = outneighbors(g, vs[u])
        if i > length(list_u)
            u += 1
            i = one(u)
            continue
        end
        if is_directed(g)
            e = MultipleEdge(vs[u], list_u[i], mul(g, vs[u], list_u[i]))
            state = (u, i + 1)
            return e, state
        else
            if list_u[i] >= vs[u]
                e = MultipleEdge(vs[u], list_u[i], mul(g, vs[u], list_u[i]))
                state = (u, i + 1)
                return e, state
            else
                i += 1
            end
        end
    end

    if n == 0 || u > n
        return nothing
    end
end

length(eit::MultipleEdgeIter) = ne(eit.g)
