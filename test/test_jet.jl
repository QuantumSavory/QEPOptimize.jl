@testitem "JET analysis" tags=[:jet] begin

using JET
using Test
using QEPOptimize

rep = report_package(QEPOptimize;
    target_modules=(QEPOptimize,),
    ignored_modules=(
        LastFrameModule(Base),
    )
)
@show rep
@test length(JET.get_reports(rep)) == 0

end
