ZXDiagram(bir::BlockIR) = convert_to_zxd(bir)
Chain(zxd::ZXDiagram) = convert_to_chain(zxd)

convert_to_gate(::Val{:X}, loc) = Gate(X, Locations(loc))
convert_to_gate(::Val{:Z}, loc) = Gate(Z, Locations(loc))
convert_to_gate(::Val{:H}, loc) = Gate(H, Locations(loc))
convert_to_gate(::Val{:S}, loc) = Gate(S, Locations(loc))
convert_to_gate(::Val{:Sdag}, loc) = Gate(AdjointOperation(S), Locations(loc))
convert_to_gate(::Val{:T}, loc) =  Gate(T, Locations(loc))
convert_to_gate(::Val{:Tdag}, loc) = Gate(AdjointOperation(S), Locations(loc))
convert_to_gate(::Val{:SWAP}, loc1, loc2) = Gate(SWAP, Locations((loc1, loc2)))
convert_to_gate(::Val{:CNOT}, loc, ctrl) = Ctrl(Gate(X, Locations(loc)), CtrlLocations(ctrl))
convert_to_gate(::Val{:CZ}, loc, ctrl) = Ctrl(Gate(Z, Locations(loc)), CtrlLocations(ctrl))
function convert_to_gate(::Val{:Rz}, loc, theta)
    if theta isa Phase
        theta = theta * π
        theta = theta.ex
    end
    return Gate(Rz(theta), Locations(loc))
end
function convert_to_gate(::Val{:Rx}, loc, theta)
    if theta isa Phase
        theta = theta * π
        theta = theta.ex
    end
    return Gate(Rx(theta), Locations(loc))
end
function convert_to_gate(::Val{:shift}, loc, theta)
    if theta isa Phase
        theta = theta * π
        theta = theta.ex
    end
    return Gate(shift(theta), Locations(loc))
end

function push_gate!(circ::Chain, gargs...) 
    push!(circ.args, convert_to_gate(gargs...))
    return circ
end
function pushfirst_gate!(circ::Chain, gargs...) 
    pushfirst!(circ.args, convert_to_gate(gargs...))
    return circ
end

function unwrap_ssa_phase(theta, ir::Core.Compiler.IRCode)
    if theta isa Core.SSAValue
        return Phase(theta, ir.stmts[theta.id][:type])
    elseif theta isa QuoteNode
        return theta.value
    elseif theta isa Core.Const
        return theta.val
    elseif theta isa Number
        return theta
    else
        error("expect SSAValue or Number")
    end
end

function canonicalize_single_location(ir::YaoHIR.Chain)
    Chain(map(canonicalize_single_location, ir.args)...)
end

function canonicalize_single_location(node::YaoHIR.Ctrl)
    if node.ctrl isa CtrlLocations && length(node.ctrl) == 1
        ctrl = CtrlLocations(node.ctrl.storage[1], node.ctrl.flags)
    else
        ctrl = node.ctrl
    end
    return YaoHIR.Ctrl(canonicalize_single_location(node.gate), node.ctrl)
end

function canonicalize_single_location(node::YaoHIR.Gate)
    length(node.locations) == 1 || return node
    return YaoHIR.Gate(node.operation, node.locations[1])
end

function convert_to_zxd(root::YaoHIR.BlockIR)
    circ = ZXDiagram(root.nqubits)

    circuit = canonicalize_single_location(root.circuit)
    for gate in YaoHIR.leaves(circuit)
        @switch gate begin
            @case Gate(&Z, loc::Locations{Int})
                push_gate!(circ, Val(:Z), plain(loc), 1//1)
            @case Gate(&X, loc::Locations{Int})
                push_gate!(circ, Val(:X), plain(loc), 1//1)
            @case Gate(&H, loc::Locations{Int})
                push_gate!(circ, Val(:H), plain(loc))
            @case Gate(&S, loc::Locations{Int})
                push_gate!(circ, Val(:Z), plain(loc), 1//2)
            @case Gate(&T, loc::Locations{Int})
                push_gate!(circ, Val(:Z), plain(loc), 1//4)
            @case Gate(shift(theta), loc::Locations{Int})
                theta = unwrap_ssa_phase(theta, root.parent)
                push_gate!(circ, Val(:Z), plain(loc), (1/π)*theta)
            @case Gate(Rx(theta), loc::Locations{Int})
                theta = unwrap_ssa_phase(theta, root.parent)
                push_gate!(circ, Val(:X), plain(loc), (1/π)*theta)
            @case Gate(Ry(theta), loc::Locations{Int})
                theta = unwrap_ssa_phase(theta, root.parent)
                push_gate!(circ, Val(:X), plain(loc),  1//2)
                push_gate!(circ, Val(:Z), plain(loc),  (1/π) * theta)
                push_gate!(circ, Val(:X), plain(loc), -1//2)
            @case Gate(Rz(theta), loc::Locations{Int})
                theta = unwrap_ssa_phase(theta, root.parent)
                push_gate!(circ, Val(:Z), plain(loc), (1/π)*theta)
            @case Gate(AdjointOperation(&S), loc::Locations{Int})
                push_gate!(circ, Val(:Z), plain(loc), 3//2)
            @case Gate(AdjointOperation(&T), loc::Locations{Int})
                push_gate!(circ, Val(:Z), plain(loc), 7//4)
            @case Ctrl(Gate(&X, loc::Locations{Int}), ctrl::CtrlLocations{Int}) # CNOT
                push_gate!(circ, Val(:CNOT), plain(loc), plain(ctrl))
            @case Ctrl(Gate(&Z, loc::Locations{Int}), ctrl::CtrlLocations{Int}) # CZ
                push_gate!(circ, Val(:CZ), plain(loc), plain(ctrl))
            @case _
                error("$gate is not supported")
        end
    end
    return circ
end

function convert_to_chain(circ::ZXDiagram{TT, P}) where {TT, P}
    spider_seq = spider_sequence(circ)
    qc = []
    for vs in spider_seq
        if length(vs) == 1
            v = vs
            q = Int(qubit_loc(circ, vs))
            push_spider_to_chain!(qc, q, phase(circ, v), spider_type(circ, v))
        elseif length(vs) == 2
            v1, v2 = vs
            q1 = Int(qubit_loc(circ, v1))
            q2 = Int(qubit_loc(circ, v2))
            push_spider_to_chain!(qc, q1, phase(circ, v1), spider_type(circ, v1))
            push_spider_to_chain!(qc, q2, phase(circ, v2), spider_type(circ, v2))
            if spider_type(circ, v1) == SpiderType.Z && spider_type(circ, v2) == SpiderType.X
                push!(qc, Ctrl(Gate(X, Locations(q2)), CtrlLocations(q1)))
            elseif spider_type(circ, v1) == SpiderType.X && spider_type(circ, v1) == SpiderType.Z
                push!(qc, Ctrl(Gate(X, Locations(q1)), CtrlLocations(q2)))
            else
                error("Spiders ($v1, $v2) should represent a CNOT")
            end
        elseif length(vs) == 3
            v1, h, v2 = vs
            spider_type(circ, h) == SpiderType.H || error("The spider $h should be a H-box")
            q1 = Int(qubit_loc(circ, v1))
            q2 = Int(qubit_loc(circ, v2))
            push_spider_to_chain!(qc, q1, phase(circ, v1), spider_type(circ, v1))
            push_spider_to_chain!(qc, q2, phase(circ, v2), spider_type(circ, v2))
            if spider_type(circ, v1) == SpiderType.Z && spider_type(circ, v2) == SpiderType.Z
                push!(qc, Ctrl(Gate(Z, Locations(q2)), CtrlLocations(q1)))
            else
                error("Spiders ($v1, $h, $v2) should represent a CZ")
            end
        else
            error("ZXDiagram's without circuit structure are not supported")
        end
    end
    return Chain(qc...)
end

function push_spider_to_chain!(qc, q, ps, st)
    if ps != 0
        if st == SpiderType.Z
            if ps == 1
                push!(qc, Gate(Z, Locations(q)))
            elseif ps == 1//2
                push!(qc, Gate(S, Locations(q)))
            elseif ps == 3//2
                push!(qc, Gate(AdjointOperation(S), Locations(q)))
            elseif ps == 1//4
                push!(qc, Gate(T, Locations(q)))
            elseif ps == 7//4
                push!(qc, Gate(AdjointOperation(T), Locations(q)))
            elseif ps != 0
                θ = ps * π
                if θ isa Phase
                    θ = θ.ex
                end
                push!(qc, Gate(shift(θ), Locations(q)))
            end
        elseif st == SpiderType.X
            if ps == 1
                push!(qc, Gate(X, Locations(q)))
            else ps != 0
                θ = ps * π
                if θ isa Phase
                    θ = θ.ex
                end
                push!(qc, Gate(Rx(θ), Locations(q)))
            end
        elseif st == SpiderType.H
            push!(qc, Gate(H, Locations(q)))
        end
    end
end