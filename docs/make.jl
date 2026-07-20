using Documenter
using QEPOptimize

makedocs(
    sitename = "QEPOptimize.jl",
    modules = [QEPOptimize],
    doctest = false,
    warnonly = :missing_docs,
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/QuantumSavory/QEPOptimize.jl.git",
    devbranch = "master",
    push_preview = true,
)
