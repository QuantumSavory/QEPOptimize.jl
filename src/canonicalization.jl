# Special thanks to Marius Minea for teaching a great class on higher-order functions! - xoth42

"""
    cleanup_two_measurements!(ops,num_pairs)

Remove unnecessary measurements (back to back on the same qubit)
"""
function cleanup_two_measurements!(ops,num_pairs)
    # This starts as [false,false...] for each pair/qubit.
    # for each op, mark this accordingly -> non-measure would set the qubits to false, and yes measure would set qubits to true. If we encounter a measure on something that already has true, it is a redundant back-to-back measure, so we delete that op. 
    last_op_on_pair_was_measure = falses(num_pairs)
    
    ops_to_delete = [] # this will collect the indexes of redundant measures to be removed

    # go through the each operation 
    for (op_index,op) in enumerate(ops) 
        measure = typeof(op) <: AbstractMeasurement # is the op a measurement?
        for qubit in affectedqubits(op)
            if measure 
                # measure found, is it redundant?
                if last_op_on_pair_was_measure[qubit]
                    # redundant measure found. add index to delete list
                    push!(ops_to_delete,op_index)
                else
                    # non redundant, but mark measured array incase there is another
                    last_op_on_pair_was_measure[qubit] = true
                end 
            else 
                # this op is not a measurement, mark measured array as false 
                last_op_on_pair_was_measure[qubit] = false 
            end
        end
    end

    # Now act on the delete list 
    return deleteat!(ops,ops_to_delete)
end

"""
    get_used_qubits(ops)

Get a Set of all qubits that are affected by the ops.
"""
function get_used_qubits(ops)
    # for each op, get its affectedqubits. (map ops to qubits)
    # affected qubits can return Int, or Vector{Int}
    # expand these results (...), making a set of all ints
    # this effectively is the set of all 'affected/used' qubits/pairs
    return reduce((set,qubits) -> push!(set, affectedqubits(qubits)...), ops; init=Set())
end

"""
    unsafe_cleanup_untargeted_pairs!(ops,num_pairs,num_purified)

UNSAFE
Check if the circuit has an unused qubit (technically a pair), if so, add a random two qubit Bell preserving gate `CNOTPerm` and coincidence measure in a random basis, ie, use it in a "probably" good way.
"""
function unsafe_cleanup_untargeted_pairs!(ops,num_pairs,num_purified)
    if num_pairs <= num_purified 
        # assumption for adding CNOTPerms.
        # if the user wants it anyway, just warn and return 
        @warn "Not enough pairs - skipping canonicalization step cleanup_untargeted_pairs!"
        return ops
    end

    # The pairs that should be used (all of them)
    pairs_to_use = Set([i for i in 1:num_pairs])
    
    # pairs that are actually used
    pairs_used = get_used_qubits(ops)

    # pairs that are unused, and others
    pairs_unused = symdiff(pairs_to_use,pairs_used)

    # use any unused pairs
    for pair in pairs_unused
        # all pairs should be in the usable pairs 
        @assert pair in pairs_to_use

        # if it is a purified pair, use this pair as a control for a new CNOT 
        if pair <= num_purified
            pushfirst!(ops,CNOTPerm(rand(1:6),rand(1:6),pair,rand(num_purified+1:num_pairs)))
        else 
            # add a CNOT, and measure to the front. New ops becomes -> [CNOT, MEAS, ...]
            # using pushfirst! so this is done in opposite order.
            # control will be a purified pair, and target is this pair
            pushfirst!(ops,rand(BellMeasure,pair))
            pushfirst!(ops,rand(CNOTPerm,rand(1:num_purified),pair))
        end
    end
    
    return ops
end

"""
    cleanup_measurements_on_top_qubits!(ops,num_purified)

Checks for and removes any measurements on purified pairs
"""
function cleanup_measurements_on_top_qubits!(ops,num_purified)
    return filter!(op-> 
        !(  # remove it if,
            typeof(op) <: AbstractMeasurement # it is a measurement and,
            && op.sidx <= num_purified # it measures a purified pair
        ), ops)
end

"""
    unsafe_cleanup_nonmeasurement_last_steps!(ops,num_pairs,num_purified)

UNSAFE
If there is a non-measurement in the last step of a non-purified pair, add random coincidence measurements.
"""
function unsafe_cleanup_nonmeasurement_last_steps!(ops,num_pairs,num_purified)
     if num_pairs <= num_purified 
        @warn "Not enough pairs - skipping canonicalization step cleanup_nonmeasurement_last_steps!"
        return ops
    end

    non_pure_last_steps = zeros(Int64,num_pairs - num_purified) # array id, plus the num_purified, is the non_pure pair number. The contents of the array 
    notDone = true 
    while notDone
        # get all of the last steps
        for (op_index,op) in enumerate(ops)
            qubits = affectedqubits(op)
            for qubit in qubits
                if qubit > num_purified
                    # this is a non purified pair, mark a new last step
                    non_pure_last_steps[qubit - num_purified] = op_index
                end
            end
        end

        # check if any of the last ops are not measurements, and ignore zeros (no ops on that qubit)
        qubits_who_need_measurements = reduce((arr, (qubit, op_index)) ->
                (op_index == 0 # filter out zeros
                || typeof(ops[op_index]) <: AbstractMeasurement) ? arr : # and filter out measurments
                push!(arr,qubit+num_purified), # otherwise, keep it and add 'num_purified' so the qubit index is correct
            enumerate(non_pure_last_steps);init=[])
        
        # all of the marked ops are non_last measurements, so add measures to them 
        for qubit in qubits_who_need_measurements
            push!(ops,rand(BellMeasure,qubit))
        end
        
        # check if done
        if length(qubits_who_need_measurements) == 0
            notDone = false
        end
    end
    return ops
end

"""
    canonicalize_cleanup!(pop::Population,num_pairs,num_purified; 
    safe=true # true to only improve circuits, false to include canonicalizations that may make some circuits worse/make experimental changes to try to improve them
    )

Apply all of the canonicalization cleanup methods to a population with tmap. Defaults to safe canonicalizations. 
"""
function canonicalize_cleanup!(pop::Population,num_pairs,num_purified; 
    safe=true # :all (include canonicalizations that may make some circuits worse/make experimental changes to try to improve) or :safe (only improve circuits)
    )
    # Function to change the ops based on the canonicalizations requested
    cleanup! = safe ?
    (indiv) -> begin # safe is true
        # Safe canonicalization
        cleanup_measurements_on_top_qubits!(indiv.ops,num_purified)
        cleanup_two_measurements!(indiv.ops,num_pairs)
    end : (indiv) -> begin # safe is false
        # Unsafe canonicalization (plus safe checks)
        cleanup_measurements_on_top_qubits!(indiv.ops,num_purified)
        cleanup_two_measurements!(indiv.ops,num_pairs)
        unsafe_cleanup_nonmeasurement_last_steps!(indiv.ops,num_pairs,num_purified)
        unsafe_cleanup_untargeted_pairs!(indiv.ops,num_pairs,num_purified)
    end 
    
    tmap(cleanup!, pop.individuals)
end
