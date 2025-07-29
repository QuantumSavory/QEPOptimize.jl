"""
    cleanup_two_measurements!(ops)

Cleanup circuits - remove unnecessary measurements
"""
function cleanup_two_measurements!(ops)
    notDone = true
    while notDone
        opRemoved = false
        # Map the list of ops to list of bools, true if that index is a measurement
        measurements = map(op -> typeof(op) <: AbstractMeasurement,ops)
        # check if there are two measurements back to back on the same register 
        for i in 1:(length(ops) -1)
            if measurements[i] && measurements[i+1] &&
                ops[i].sidx == ops[i+1].sidx
                # back to back measurement on the same pair 
                # Remove the last one 
                deleteat!(ops,i+1)
                # now restart the outer loop
                opRemoved = true
                break
            end
        end

        # finished a scan, did we find something?
        if !opRemoved
            # did not find anything, we are done
            notDone = false
        end
    end
    return ops
end
