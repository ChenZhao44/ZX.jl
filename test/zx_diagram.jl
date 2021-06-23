g = Multigraph([0 1 0; 1 0 1; 0 1 0])
ps = [Rational(0) for i = 1:3]
v_t = [SpiderType.X, SpiderType.Z, SpiderType.X]
zxd = ZXDiagram(g, v_t, ps)
zxd2 = ZXDiagram(g, Dict(zip(1:3,v_t)), Dict(zip(1:3,ps)))
@test zxd.mg == zxd2.mg && zxd.st == zxd2.st && zxd.ps == zxd2.ps

zxd2 = copy(zxd)
@test zxd.st == zxd2.st && zxd.ps == zxd2.ps
@test ZXCalculus.spider_type(zxd, 1) == SpiderType.X
@test nv(zxd) == 3 && ne(zxd) == 2

@test rem_edge!(zxd, 2, 3)
@test outneighbors(zxd, 2) == inneighbors(zxd, 2)

ZXCalculus.add_spider!(zxd, SpiderType.H, 0//1, [2, 3])
ZXCalculus.insert_spider!(zxd, 2, 4, SpiderType.H)
@test nv(zxd) == 5 && ne(zxd) == 4

zxd3 = ZXDiagram(3)
ZXCalculus.insert_spider!(zxd3, 1, 2, SpiderType.H)
pushfirst_gate!(zxd3, Val{:SWAP}(), [1, 2])
push_gate!(zxd3, Val{:SWAP}(), [2, 3])
@test ZXCalculus.qubit_loc(zxd3, 1) == ZXCalculus.qubit_loc(zxd3, 2)

@testset "float to rational" begin
    @test ZXCalculus.continued_fraction(2.41, 10) === 241//100
    @test ZXCalculus.continued_fraction(1.3, 10) === 13//10
    @test ZXCalculus.continued_fraction(0, 10) === 0//1
    @test ZXCalculus.continued_fraction(-0.5, 10) === -1//2
    zxd = ZXDiagram(4)
    push_gate!(zxd, Val(:X), 3, 0.5)
    @test zxd.ps[9] == 1//2
    push_gate!(zxd, Val(:X), 3, -0.5)
    @test zxd.ps[10] == 3//2
    push_gate!(zxd, Val(:Z), 3, 0)
    @test zxd.ps[11] == 0//1
    @test_warn "" push_gate!(zxd, Val(:Z), 3, sqrt(2))
    @test_throws MethodError push_gate!(zxd, Val(:Z), 3, sqrt(2); autoconvert=false)
    @test ZXCalculus.safe_convert(Rational{Int64}, 1.2) == 6//5 && ZXCalculus.safe_convert(Rational{Int64}, 1//2) == 1//2
end

zxd4 = ZXDiagram(2)
ZXCalculus.add_global_phase!(zxd4, ZXCalculus.Phase(1//2))
ZXCalculus.add_power!(zxd4, 2)
@test ZXCalculus.scalar(zxd4) == ZXCalculus.Scalar(2, 1//2)
pushfirst_gate!(zxd4, Val(:X), 1)
pushfirst_gate!(zxd4, Val(:H), 1)
pushfirst_gate!(zxd4, Val(:CNOT), 2, 1)
pushfirst_gate!(zxd4, Val(:CZ), 1, 2)
@test indegree(zxd4, 5) == outdegree(zxd4, 5) == degree(zxd4, 5)