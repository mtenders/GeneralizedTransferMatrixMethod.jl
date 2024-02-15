cd(@__DIR__)
push!(LOAD_PATH,"../src/")

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Documenter
using GeneralizedTransferMatrixMethod

makedocs(
    sitename="GeneralizedTransferMatrixMethod.jl",
    authors = "Michael T. Enders",
    pages = [
        "Introduction" => "index.md"
        "Manual" => Any[
            "gettingstarted.md",
            "unitful.md",
            # "parallelization.md",
            # "examples.md"
        ]
        "Library" => Any[
            "library-public.md",
            "library-internal.md"
        ]
    ],
)

deploydocs(
    repo = "github.com/mtenders/GeneralizedTransferMatrixMethod.jl",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"]
)
