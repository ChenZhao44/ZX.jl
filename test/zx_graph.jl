g = Multigraph(6)
add_edge!(g, 1, 3)
add_edge!(g, 2, 4)
add_edge!(g, 3, 4)
add_edge!(g, 3, 5)
add_edge!(g, 4, 6)
ps = [0//1 for i = 1:6]
v_t = [SpiderType.In, SpiderType.In, SpiderType.X, SpiderType.Z, SpiderType.Out, SpiderType.Out]
zxd = ZXDiagram(g, v_t, ps)
zxg1 = ZXGraph(zxd)
@test outneighbors(zxg1, 1) == inneighbors(zxg1, 1)
@test !ZXCalculus.is_hadamard(zxg1, 2, 4) && !ZXCalculus.is_hadamard(zxg1, 4, 6)
@test add_edge!(zxg1, 1, 1)
@test !add_edge!(zxg1, 2, 4)
@test !add_edge!(zxg1, 7, 8)
@test sum([ZXCalculus.is_hadamard(zxg1, src(e), dst(e)) for e in edges(zxg1.mg)]) == 3
replace!(Rule{:b}(), zxd)
zxg2 = ZXGraph(zxd)
@test !ZXCalculus.is_hadamard(zxg2, 5, 8) && !ZXCalculus.is_hadamard(zxg2, 1, 7)

zxd = ZXDiagram(2)
push_gate!(zxd, Val(:H), 1)
push_gate!(zxd, Val(:CNOT), 2, 1)
zxg = ZXGraph(zxd)

zxg3 = ZXGraph(ZXDiagram(3))
ZXCalculus.add_global_phase!(zxg3, ZXCalculus.Phase(1//4))
ZXCalculus.add_power!(zxg3, 3)
@test ZXCalculus.scalar(zxg3) == Scalar(3, 1//4)
@test degree(zxg3, 1) == indegree(zxg3, 1) == outdegree(zxg3, 1)
@test ZXCalculus.qubit_loc(zxg3, 1) == ZXCalculus.qubit_loc(zxg3, 2)
@test ZXCalculus.column_loc(zxg3, 1) == 1//1
@test ZXCalculus.column_loc(zxg3, 2) == 3//1