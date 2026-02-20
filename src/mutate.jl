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

"make a 'child' individual from meshing together two individuals' operations"
function make_child(mother_ops::Vector{Any}, father_ops::Vector{Any},max_ops::Int64)
    # one/both empty ops edge case
    mother_empty = length(mother_ops) == 0
    father_empty = length(father_ops) == 0
    mother_empty && !father_empty && return Individual(:child, copy(father_ops))
    !mother_empty && father_empty && return Individual(:child, copy(mother_ops))
    mother_empty && father_empty && return Individual(:child, [])

    # Child algorithm:
    # choose a location to 'split' the circuits of both the mother and father. The mother and father have their own splits. All ops up to the split will be taken from the mother, and all ops from the end-minus-split to the end will be taken from the father. 

    # to make sure that we do not pass the max ops limit, the max split will be half of the max ops. This is because the amount of ops in the end will be mother_split + father_split
    mother_split = rand(1:min(convert(Int,floor(max_ops/2)),length(mother_ops)))
    father_split = rand(1:min(convert(Int,floor(max_ops/2)),length(father_ops)))
    
    # mark the history, and concat the ops
    child = Individual(:child,[mother_ops[1:mother_split]; father_ops[end-father_split+1:end]])
    @assert (mother_split + father_split) == length(child.ops)
    return child
end

"make a new individual with two operations randomly swapped"
function swap_op(ops::Vector{Any})
    @assert length(ops) >= 2 

    # select locations
    l = length(ops)
    i = rand(1:l)
    j = i 
    while (j==i) j = rand(1:l); end # select until i != j  

    # swap the ops
    new_ops = copy(ops) 
    new_ops[j],new_ops[i] = new_ops[i],new_ops[j]
    return Individual(:swap,new_ops)
end
