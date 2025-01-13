"perform a mutation"
function mutate end

function mutate(op) # a default no-op "mutation" in case one is not implemented
    return op
end

function mutate(indiv::Individual)
    if length(indiv.ops) == 0
        return indiv
    end
    mutated_ops = mutate.(indiv.ops)
    return Individual(:opsmutate, mutated_ops)
end

function mutate(gate::BellMeasure)
    return rand(BellMeasure, gate.sidx) # TODO (low priority) this `rand` is a "pun"; should be changed to use a keyword argument to specify affected qubit, but that is a breaking change in BPGates.jl
end

function mutate(gate::CNOTPerm)
    return rand(CNOTPerm, gate.idx1, gate.idx2) # TODO (low priority) this `rand` is a "pun"; should be changed to use a keyword argument to specify affected qubit, but that is a breaking change in BPGates.jl
end


is_droppable(::Any) = false
is_droppable(::CNOTPerm) = true
is_droppable(::BellMeasure) = true

"make a new individual with randomly deleted operations"
function drop_op(indiv::Individual)
    # Filter the indices of operations that can be dropped
    drop_indices = [i for (i,op) in pairs(indiv.ops) if is_droppable(op)]

    if  isempty(drop_indices)
        # If there are no droppable operations, return the individual as is
        return copy(indiv)
    else
        # Randomly select and delete one of the droppable operations
        new_ops = copy(indiv.ops)
        deleteat!(new_ops, rand(drop_indices))
        return Individual(:drop, new_ops)
    end
end


"make a new individual with a randomly gained operation"
function gain_op(indiv::Individual; valid_pairs)
    new_ops = copy(indiv.ops)

    op = rand_op(valid_pairs)

    position = rand(1:length(new_ops)+1)

    insert!(new_ops, position, op)

    return Individual(:gain, new_ops)
end

function rand_op(valid_pairs)
    # weighted randomly select a CNOTPerm or a measurement # TODO (low priority) make the selection and weights configurable
    op = if rand() < 0.7 && length(valid_pairs) >= 1
        i1, i2 = randperm(length(valid_pairs))[1:2]
        pair1, pair2 = valid_pairs[i1], valid_pairs[i2]
        random_gate = rand(CNOTPerm, pair1, pair2)
    else
        pair = rand(valid_pairs)
        rand(BellMeasure, pair)
    end
    return op
end
