# Unitful.jl integration

```@setup unitful
using Plots

default(
    lw=2, 
    label=:none,
    framestyle=:box,
    grid=false
)
```

GeneralizedTransferMatrixMethod.jl automatically loads and reexports
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl). This way all input
parameters can be specified using the units supported by Unitful.jl


```@example unitful
using LinearAlgebra

using GeneralizedTransferMatrixMethod
# Import most common units 
using Unitful: °, nm, μm, mm, m
```

We simulate a slab of a uniaxial crystal in air as an illustrative example
```@example unitful
Air = Layer()

nₒ(λ) = 3
nₑ(λ) = 1.5

@permittivity "Mat" λ -> Diagonal([nₒ(λ), nₑ(λ), nₒ(λ)].^2)
nothing # hide
```

We can define all input parameters using Unitful.jl quantities
```@example unitful
α = 30°
ζ = sin(α)
ϕ = 45°
d = 183nm

λ = (500:0.1:600)nm
nothing # hide
```

and pass them to setup the [`LayeredStructure`](@ref LayeredStructure)
```@example unitful
Stack = LayeredStructure(
    superstrate = Air,
    layers = [Mat(d = d, ϕ = ϕ)],
    substrate = Air
)
nothing # hide
```

and use them to calculate its properties
```@example unitful
Properties = [calculate_structure_properties(ζ, λᵢ, Stack) for λᵢ ∈ λ]

using Unzip

Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ = unzip(reflection.(Properties))

T(prop) = abs2.(transmission_coeffs(prop))
Tₚₚ, Tₛₛ, Tₚₛ, Tₛₚ = unzip(T.(Properties))
nothing # hide
```

Units get automatically added to the axis labels in Plots.jl
```@example unitful
using Plots

p1 = plot(
    λ, [Rₚₚ Rₛₛ Rₚₛ Rₛₚ], 
    label = ["Rₚₚ" "Rₛₛ" "Rₚₛ" "Rₛₚ"],
    xlabel = "Wavelength",
    ylabel = "Reflection",
    ls = [:solid :dash],
)
p2 = plot(
    λ, [Tₚₚ Tₛₛ Tₚₛ Tₛₚ], 
    label = ["Tₚₚ" "Tₛₛ" "Tₚₛ" "Tₛₚ"],
    xlabel = "Wavelength",
    ylabel = "Transmission",
    ls = [:solid :dash],
)
plot(p1,p2)
```
