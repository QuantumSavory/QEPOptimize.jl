using QEPOptimize
using TestItemRunner

# filter for the test
testfilter = ti -> begin
  exclude = Symbol[]
  # From QuantumSavory.jl
  if get(ENV,"JET_TEST","")!="true"
    push!(exclude, :jet)
  end

  return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter

