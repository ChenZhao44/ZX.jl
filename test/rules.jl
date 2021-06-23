g = Multigraph([0 2 0; 2 0 1; 0 1 0])
collect(edges(g))
ps = [i//4 for i = 1:3]
v_t = [SpiderType.Z, SpiderType.Z, SpiderType.X]
zxd = ZXDiagram(g, v_t, ps)
matches = match(Rule{:f}(), zxd)
rewrite!(Rule{:f}(), zxd, matches)
@test sort!(spiders(zxd)) == [1, 3]
@test phase(zxd, 1) == phase(zxd, 3) == 3//4

g = Multigraph(path_graph(5))
add_edge!(g, 1, 2)
ps = [1, 0//1, 0, 0, 1]
v_t = [SpiderType.X, SpiderType.X, SpiderType.Z, SpiderType.Z, SpiderType.Z]
zxd = ZXDiagram(g, v_t, ps)
matches = match(Rule{:i1}(), zxd)
rewrite!(Rule{:i1}(), zxd, matches)
@test nv(zxd) == 3 && ne(zxd, count_mul = true) == 3 && ne(zxd) == 2

g = Multigraph([0 2 0; 2 0 1; 0 1 0])
ps = [i//4 for i = 1:3]
v_t = [SpiderType.X, SpiderType.X, SpiderType.Z]
zxd = ZXDiagram(g, v_t, ps)
matches = match(Rule{:h}(), zxd)
rewrite!(Rule{:h}(), zxd, matches)
@test nv(zxd) == 8 && ne(zxd) == 8

matches = match(Rule{:i2}(), zxd)
rewrite!(Rule{:i2}(), zxd, matches)
@test nv(zxd) == 4 && ne(zxd, count_mul = true) == 4 && ne(zxd) == 3

g = Multigraph(6)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 3, 5)
add_edge!(g, 3, 6)
ps = [0, 1, 1//2, 0, 0, 0]
v_t = [SpiderType.In, SpiderType.X, SpiderType.Z, SpiderType.Out, SpiderType.Out, SpiderType.Out]
zxd = ZXDiagram(g, v_t, ps)
matches = match(Rule{:pi}(), zxd)
rewrite!(Rule{:pi}(), zxd, matches)
@test nv(zxd) == 8 && ne(zxd) == 7
@test zxd.scalar == Scalar(0, 1//2)

g = Multigraph([0 2 0; 2 0 1; 0 1 0])
ps = [1, 1//2, 0]
v_t = [SpiderType.X, SpiderType.Z, SpiderType.In]
zxd = ZXDiagram(g, v_t, ps)
matches = match(Rule{:pi}(), zxd)
rewrite!(Rule{:pi}(), zxd, matches)
@test nv(zxd) == 4 && ne(zxd) == 3 && ne(zxd, count_mul = true) == 4
@test zxd.scalar == Scalar(0, 1//2)

g = Multigraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3, 2)
add_edge!(g, 2, 4)
add_edge!(g, 2, 5)
ps = [0, 1//2, 0, 0, 0]
v_t = [SpiderType.X, SpiderType.Z, SpiderType.Out, SpiderType.Out, SpiderType.Out]
zxd = ZXDiagram(g, v_t, ps)
matches = match(Rule{:c}(), zxd)
rewrite!(Rule{:c}(), zxd, matches)
@test nv(zxd) == 7 && ne(zxd) == 4
@test zxd.scalar == Scalar(-3, 0//1)

g = Multigraph(6)
add_edge!(g, 1, 3)
add_edge!(g, 2, 4)
add_edge!(g, 3, 4)
add_edge!(g, 3, 5)
add_edge!(g, 4, 6)
ps = [0//1 for i = 1:6]
v_t = [SpiderType.In, SpiderType.In, SpiderType.X, SpiderType.Z, SpiderType.Out, SpiderType.Out]
layout = ZXCalculus.ZXLayout(2, Dict(zip(1:6, [1//1, 2, 1, 2, 1, 2])), Dict(zip(1:6, [1//1, 1, 2, 2, 3, 3])))
zxd = ZXDiagram(g, v_t, ps, layout)
matches = match(Rule{:b}(), zxd)
rewrite!(Rule{:b}(), zxd, matches)
@test nv(zxd) == 8 && ne(zxd) == 8
@test zxd.scalar == Scalar(1, 0//1)

g = Multigraph(9)
for e in [[2,6],[3,7],[4,8],[5,9]]
    add_edge!(g, e[1], e[2])
end
ps = [1//2, 0, 1//4, 1//2, 3//4, 0, 0, 0, 0]
st = [SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.In, SpiderType.In, SpiderType.Out, SpiderType.Out]
zxg = ZXGraph(ZXDiagram(g, st, ps))
for e in [[1,2],[1,3],[1,4],[1,5],[2,3]]
    add_edge!(zxg, e[1], e[2])
end
replace!(Rule{:lc}(), zxg)
@test !has_edge(zxg, 2, 3) && ne(zxg) == 9
@test phase(zxg, 2) == 3//2 && phase(zxg, 3) == 7//4 && phase(zxg, 4) == 0//1 && phase(zxg, 5) == 1//4

g = Multigraph(14)
for e in [[3,9],[4,10],[5,11],[6,12],[7,13],[8,14]]
    add_edge!(g, e[1], e[2])
end
ps = [1//1, 0, 1//4, 1//2, 3//4, 1, 5//4, 3//2, 0, 0, 0, 0, 0, 0]
st = [SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.In, SpiderType.Out, SpiderType.In, SpiderType.Out, SpiderType.In, SpiderType.Out]
zxg = ZXGraph(ZXDiagram(g, st, ps))
for e in [[1,2],[1,3],[1,4],[1,5],[1,6],[2,5],[2,6],[2,7],[2,8]]
    add_edge!(zxg, e[1], e[2])
end

replace!(Rule{:p1}(), zxg)
@test !has_edge(zxg, 3, 4) && !has_edge(zxg, 5, 6) && !has_edge(zxg, 7, 8)
@test nv(zxg) == 12 && ne(zxg) == 18
@test phase(zxg, 3) == 1//4 && phase(zxg, 4) == 1//2 && phase(zxg, 5) == 3//4 && phase(zxg, 6) == 1//1 && phase(zxg, 7) == 1//4 && phase(zxg, 8) == 1//2

g = Multigraph(6)
for e in [[2,6]]
    add_edge!(g, e[1], e[2])
end
ps = [1//1, 1//4, 1//2, 3//4, 1, 0]
st = [SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.In]
zxg = ZXGraph(ZXDiagram(g, st, ps))
for e in [[1,2],[2,3],[1,4],[1,5]]
    add_edge!(zxg, e[1], e[2])
end

@test length(match(Rule{:p1}(), zxg)) == 1
replace!(Rule{:pab}(), zxg)
@test nv(zxg) == 6 && ne(zxg) == 6

g = Multigraph(14)
for e in [[3,9],[4,10],[5,11],[6,12],[7,13],[8,14]]
    add_edge!(g, e[1], e[2])
end
ps = [1//1, 1//4, 1//4, 1//2, 3//4, 1, 5//4, 3//2, 0, 0, 0, 0, 0, 0]
st = [SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z,
    SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.In, SpiderType.Out,
    SpiderType.In, SpiderType.Out, SpiderType.In, SpiderType.Out]
zxg = ZXGraph(ZXDiagram(g, st, ps))
for e in [[1,2],[1,3],[1,4],[1,5],[1,6],[2,5],[2,6],[2,7],[2,8]]
    add_edge!(zxg, e[1], e[2])
end
match(Rule{:p2}(), zxg)
replace!(Rule{:p2}(), zxg)
@test zxg.phase_ids[15] == (2, -1)

g = Multigraph(15)
for e in [[3,9],[4,10],[5,11],[6,12],[7,13],[8,14],[2,15]]
    add_edge!(g, e[1], e[2])
end
ps = [1//1, 1//4, 1//2, 1//2, 3//2, 1, 1//2, 3//2, 0, 0, 0, 0, 0, 0, 0, 0]
st = [SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.Z,
    SpiderType.Z, SpiderType.Z, SpiderType.Z, SpiderType.In, SpiderType.Out,
    SpiderType.In, SpiderType.Out, SpiderType.In, SpiderType.Out, SpiderType.Out]
zxg = ZXGraph(ZXDiagram(g, st, ps))
for e in [[1,2],[1,3],[1,4],[1,5],[1,6],[2,5],[2,6],[2,7],[2,8]]
    add_edge!(zxg, e[1], e[2])
end
replace!(Rule{:p3}(), zxg)

@test nv(zxg) == 16 && ne(zxg) == 28
@test ZXCalculus.is_hadamard(zxg, 2, 15) && ZXCalculus.is_hadamard(zxg, 1, 16)
