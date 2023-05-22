# Getting Started

This section offers a concise and user-friendly guide to effectively utilize
[GeneralizedTransferMatrixMethod.jl](https://github.com/mtenders/GeneralizedTransferMatrixMethod.jl). 

## Installation

The latest version of GeneralizedTransferMatrixMethod.jl can be installed via
the built-in package manager. In the Julia REPL, press `]` to access the package
manager and run
```
pkg> add GeneralizedTransferMatrixMethod
```
Updates can be installed similarly using the `up` command
```
pkg> up GeneralizedTransferMatrixMethod
```
!!! warn "Julia version"
    It is recommended to use a recent Julia version. This package was tested on
    Julia v1.7 and onwards, but does not provide any backwards compatibility.

## Tutorial

This tutorial goes through some simple examples to explain the key features of 
GeneralizedTransferMatrixMethod.jl.

!!! note "Unitful.jl integration"
    In this tutorial we all input parameters are given in SI units. This package
    supports the use of arbitrary units using the
    [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) packages. For more
    information see [Unitful.jl integration](@ref).

### Glass--Air interface

```@setup tutorial1
using GeneralizedTransferMatrixMethod
using Plots

default(
    lw=2, 
    label=:none,
    framestyle=:box,
    grid=false
)
```

We should now be able to load the package

```julia
using GeneralizedTransferMatrixMethod
```

We first have to define the materials we want to use. While there are a few
materials included in the package ([`Au`](@ref Au), [`Ag`](@ref Ag),
[`SiC`](@ref SiC), [`MoO₃`](@ref MoO₃), ...), we start by defining the
permittivities from scratch. 

For this we use the [`@permittivity`](@ref @permittivity) macro. This macro
expects the name of the material we want to define and a function that takes the
wavelength `λ` in meters and returns the full permittivity tensor.

A simple example would be to define Glass with a constant permittivity of 1.5
```@example tutorial1
using LinearAlgebra

@permittivity "Glass" λ -> Diagonal(ones(3)) * 1.5
nothing #hide
```
We first load the `LinearAlgebra` standard library to have access to the
function `Diagonal`, so we don't have to write out the ``3 \times 3`` diagonal
matrix ourselves. In this example we passed a anonymous function (`->` syntax)to
the permittivity macro, but we could have also defined the function beforehand
```julia
f(λ) = Diagonal(ones(3)) * 1.5

@permittivity "Glass" f
```

The macro has defined a function called "Glass," which returns a structure of
type [`Layer`](@ref Layer). Additionally, it has created a function `ϵ_Glass`,
which is identical to the input function we provided.

```@example tutorial1
Glass()
```

To calculate the reflection of a Glass--Air interface, we must define an Air
layer. To accomplish this, we can simply invoke the [`Layer`](@ref Layer)
function without providing any arguments. This will create a Layer with a
constant permittivity of 1
```@example tutorial1
Air = Layer()
```

With these two layers, we can now build up a [`LayeredStructure`](@ref
LayeredStructure) 
```@example tutorial1
Interface = LayeredStructure(
    superstrate = Glass(),
    substrate = Air
)
```

We can now calculate the transfer matrix of our interface for a given wavelength
`λ` and angle of incidence `α`
```@example tutorial1
λ = 1.55e-6 # [m]
α = deg2rad(10) # [rad]
ζ = sin(α) * sqrt(ϵ_Glass(λ)[1,1])

Properties = calculate_structure_properties(ζ, λ, Interface)
```
Instead of the angle of incidence `α` the function
[`calculate_structure_properties`](@ref calculate_structure_properties) actually
takes the reduced in-plane wave vector `ζ` as an input. This structure of type
[`StructureProperties`](@ref StructureProperties) contains the full transfer
matrix `Γ` (denoted as ``\Gamma^*`` in reference
[[3](https://doi.org/10.1364/JOSA.62.000502)]). 
```@example tutorial1
Properties.Γ
```
It also contains all the parameters calculated along the way for each individual
layer
```@example tutorial1
Properties.superstrate
```

We can now calculate the [`reflection`](@ref reflection) and
[`transmission`](@ref transmission) from these properties
```@example tutorial1
Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ = reflection(Properties)
```
```@example tutorial1
Tₚ, Tₛ = transmission(ζ, Properties)
```
!!! danger "Transmission modes"
    Since it is not generally possible to separate the modes in birefringent
    media, the [`transmission`](@ref transmission) function calculates the
    transmission only considering the input polarization (for more details see
    reference [[5](http://doi.org/10.1103/PhysRevB.101.165425)]). In cases where
    it is possible, one can use the [`transmission_coeffs`](@ref
    transmission_coeffs) function to calculate the transmission.

We can now calculate the angular dependence of the reflection of our interface
by replacing the incident angle `α` by a list of angles
```@example tutorial1
λ = 1.55e-6 # [m]
α = deg2rad.(0:0.1:89) # [rad]
ζ = sin.(α) * sqrt(ϵ_Glass(λ)[1,1])

Properties = calculate_structure_properties.(ζ, Ref(λ), Ref(Interface))
R = reflection.(Properties)
nothing #hide
```
!!! note "Broadcasting"
    To make sure we only broadcast over the reduced in-plane wave vectors `ζ` we
    use the `Ref` command to specify which quantities should not be used for
    broadcasting. 

We now have a list of tupels `R`, where each element has the shape `(Rₚₚ, Rₛₛ,
Rₚₛ, Rₛₚ)`. If we want to split them into separate lists, we can use the
[Unzip.jl](https://github.com/bramtayl/Unzip.jl) package
```@example tutorial1
using Unzip

Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ = unzip(R)
nothing #hide
```

Finally to visualize our results, we use
[Plots.jl](https://github.com/JuliaPlots/Plots.jl/)
```@example tutorial1
using Plots

plot(
    rad2deg.(α), [Rₚₚ Rₛₛ], 
    label = ["Rₚₚ" "Rₛₛ"],
    xlabel = "Angle of incidence (°)", ylabel = "Reflection"
)
```

### Half-wave plate with anti-reflection coating
