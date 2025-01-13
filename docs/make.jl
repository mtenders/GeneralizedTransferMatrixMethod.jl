cd(@__DIR__)
push!(LOAD_PATH,"../src/")

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Documenter
using GeneralizedTransferMatrixMethod


DocMeta.setdocmeta!(GeneralizedTransferMatrixMethod, :DocTestSetup, :(using
    GeneralizedTransferMatrixMethod); recursive=true)
makedocs(
    sitename="GeneralizedTransferMatrixMethod.jl",
    authors = "Michael T. Enders",
    doctest = true,
    modules = [GeneralizedTransferMatrixMethod],
    pages = [
        "Introduction" => "index.md",
        "Manual" => Any[
             "gettingstarted.md",
             "unitful.md"
        ],
        "Examples" => Any[
            "sic-sphp.md",
            "photonic-crystal.md",
            "nonreciprocity.md"
        ],
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
