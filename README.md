<div align="center">
    
![header](./docs/src/assets/banner.svg)

---

[![DOCS](https://img.shields.io/badge/docs-GeneralizedTransferMatrixMethod.jl-blue?style=flat-square)](https://mtenders.github.io/GeneralizedTransferMatrixMethod.jl/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.7974657-blue?style=flat-square)](https://doi.org/10.5281/zenodo.7974657)

</div>

This package implements the transfer matrix method  as presented in
[[1](#References)] for arbitrary relative permittivities ``ϵ``,
*including non-reciprocal cases*. Arbitrary relative permeabilities ``μ`` and
optical tensors ``ξ`` and ``χ`` are also implemented but not yet tested.

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add GeneralizedTransferMatrixMethod
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("GeneralizedTransferMatrixMethod")
```


## Example

This example of an glass–air interface is taken from the tutorial presented in the documentation

```julia
using GeneralizedTransferMatrixMethod
using LinearAlgebra
using Unitful: °, nm, μm, mm, m

# Define Materials
@permittivity "Glass" λ -> Diagonal(ones(3)) * 1.5^2
Air = Layer()

Interface = LayeredStructure(
    superstrate = Glass(),
    substrate = Air
)

# Calculate properties
λ = 1.55μm
α = 10°

Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ = calculate_reflection(λ, α, Interface)
Tₚₚ, Tₛₛ, Tₚₛ, Tₛₚ = calculate_transmission(λ, α, Interface)
```

## More Examples
You can find more examples of relevant use causes in the [documentation](https://mtenders.github.io/GeneralizedTransferMatrixMethod.jl/). These examples are also available as [Pluto](https://github.com/fonsp/Pluto.jl) notebooks:

#### 👉 [Surface phonon polariton in 6H-SiC](https://mtenders.github.io/transfer-matrix-examples/surface_phonon_polariton.html)
#### 👉 [Distributed Bragg reflector](https://mtenders.github.io/transfer-matrix-examples/quarter_wave_stack.html)
#### 👉 More to come ...

If you want to run these notebooks on your own computer press "Edit or run this notebook" in the top right corner and follow the instructions.

## References
1. [Mackay, T. G. & Lakhtakia, A. The Transfer-Matrix Method in Electromagnetics and Optics. vol. 1 (2020).](https://doi.org/10.1007/978-3-031-02022-3)
