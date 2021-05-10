using Test
using YaoHIR
using YaoLocations
using CompilerPluginTools
using ZXCalculus
using ZXCalculus: Phase, BlockIR
using YaoHIR: Chain, Gate, Ctrl, shift, Rz

ir = @make_ircode begin
    Expr(:call, :+, 1, 1)::Int
    Expr(:call, :+, 3, 3)::Int
end

circ = Chain(
    Gate(shift(1.0), Locations(1)),
    Gate(Rz(Core.SSAValue(2)), Locations(1)),
)

bir = BlockIR(ir, 4, circ)

zxd = convert_to_zxd(bir)
qc_tl = convert_to_block_ir(phase_teleportation(zxd))
@test length(qc_tl) == 1
