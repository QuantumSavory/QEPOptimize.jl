
# TODO
@testitem "to_stabilizer runs" setup=[setup] begin
    using QEPOptimize:to_stabilizer

    output = to_stabilizer(setup.test_ops_small,4)

    to_stabilizer(setup.test_ops_small,4;show_steps=true)
end
