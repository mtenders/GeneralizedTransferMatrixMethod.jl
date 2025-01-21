# Unitful.jl integration


```@setup unitful
using Plots

default(
    lw=3, 
    label=:none,
    framestyle=:box,
    grid=false,
    size = (800,300),
    bottom_margin=5Plots.mm,
    left_margin=5Plots.mm,
    right_margin=5Plots.mm
)
```

GeneralizedTransferMatrixMethod.jl has an extension for
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl), which is loaded if
both packages are used together. This way all input parameters can be specified
using the units supported by Unitful.jl.


```@example unitful
using LinearAlgebra

using GeneralizedTransferMatrixMethod
# Import most common units 
using Unitful: °, nm, μm, mm, m
```

!!! tip "Wavenumber, Frequency, Energy"
    The use wavenumber units (e.g. `u"cm^-1"`), frequency units (e.g. `u"Hz"`)
    and energy units (e.g. `u"eV"`) is also supported out of the box.

We simulate a slab of a uniaxial crystal in air as an illustrative example
```@example unitful
Air = Layer()

nₒ(λ) = 3
nₑ(λ) = 1.5

@permittivity "Mat" λ -> Diagonal([nₒ(λ), nₑ(λ), nₒ(λ)].^2)
```

We can define all input parameters using Unitful.jl quantities
```@example unitful
α = 30°
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
```

and use them to calculate its properties
```@example unitful

using Unzip

R = calculate_reflection.(λ, α, Stack)
Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ = unzip(R)

T = calculate_transmission.(λ, α, Stack)
Tₚₚ, Tₛₛ, Tₚₛ, Tₛₚ = unzip(T)
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
savefig("unitful-plot.svg"); nothing # hide
```

![](unitful-plot.svg)
