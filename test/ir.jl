using Test
using ZXCalculus
using YaoHIR, YaoLocations
using YaoHIR.IntrinsicOperation
using CompilerPluginTools

chain = Chain()
push_gate!(chain, Val(:Sdag), 1)
push_gate!(chain, Val(:Z), 1)
push_gate!(chain, Val(:H), 1)
push_gate!(chain, Val(:S), 1)
push_gate!(chain, Val(:S), 2)
push_gate!(chain, Val(:H), 4)
push_gate!(chain, Val(:CNOT), 3, 2)
push_gate!(chain, Val(:CZ), 4, 1)
push_gate!(chain, Val(:H), 2)
push_gate!(chain, Val(:T), 2)
push_gate!(chain, Val(:CNOT), 3, 2)
push_gate!(chain, Val(:Tdag), 2)
push_gate!(chain, Val(:CNOT), 1, 4)
push_gate!(chain, Val(:H), 1)
push_gate!(chain, Val(:T), 2)
push_gate!(chain, Val(:S), 3)
push_gate!(chain, Val(:H), 4)
push_gate!(chain, Val(:T), 1)
push_gate!(chain, Val(:T), 2)
push_gate!(chain, Val(:T), 3)
push_gate!(chain, Val(:X), 3)
push_gate!(chain, Val(:H), 2)
push_gate!(chain, Val(:H), 3)
push_gate!(chain, Val(:Sdag), 4)
push_gate!(chain, Val(:S), 3)
push_gate!(chain, Val(:X), 4)
push_gate!(chain, Val(:CNOT), 3, 2)
push_gate!(chain, Val(:H), 1)
push_gate!(chain, Val(:shift), 4, ZXCalculus.Phase(1//2))
push_gate!(chain, Val(:Rx), 4, ZXCalculus.Phase(1//1))
push_gate!(chain, Val(:Rx), 3, ZXCalculus.Phase(1//4))
push_gate!(chain, Val(:Rx), 2, ZXCalculus.Phase(1//4))
push_gate!(chain, Val(:S), 3)

ir = @make_ircode begin
end
bir = BlockIR(ir, 4, chain)
zxd = convert_to_zxd(bir)
convert_to_chain(zxd)
pt_zxd = phase_teleportation(zxd)
@test tcount(pt_zxd) <= tcount(zxd)
pt_chain = convert_to_chain(pt_zxd)
@test length(pt_chain) <= length(chain)

zxg = clifford_simplification(zxd)
cl_chain = circuit_extraction(zxg)

zxg = full_reduction(zxd)
fl_chain = circuit_extraction(zxg)
ZXCalculus.generate_layout!(zxg)
@test ZXCalculus.qubit_loc(zxg, 40) == 0//1
ZXCalculus.spider_sequence(zxg)

pt_bir = phase_teleportation(bir)
cl_bir = clifford_simplification(bir)
fl_bir = full_reduction(bir)

@test length(pt_chain) == length(pt_bir.circuit)
@test length(cl_chain) == length(cl_bir.circuit)
@test length(fl_chain) == length(fl_bir.circuit)

@testset "issue#80" begin
    ir = @make_ircode begin
    end

    circuit = Chain(Gate(X, Locations(1)), Gate(X, Locations(1)))
    bir = BlockIR(ir, 1, circuit)
    bir = clifford_simplification(bir)
    bir = clifford_simplification(bir)
    @test bir.circuit == Chain(Gate(H, Locations((1, ))), Gate(H, Locations((1, ))))
end
