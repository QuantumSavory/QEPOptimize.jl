@testitem "Aqua analysis" tags=[:aqua] begin

using Aqua, QEPOptimize

Aqua.test_all(QEPOptimize)

end
