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
    Julia v1.11 and onwards, but does not provide any backwards compatibility.

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
    lw=3, 
    label=:none,
    framestyle=:box,
    grid=false,
    size = (1200,600),
    bottom_margin=5Plots.mm,
    left_margin=5Plots.mm,
    right_margin=5Plots.mm
)
```

We should now be able to load the package:
```julia
using GeneralizedTransferMatrixMethod
```

#### Setting up the structure

We first have to define the materials we want to use. While there are a few
materials included in the package ([`Au`](@ref Au), [`Ag`](@ref Ag),
[`SiC`](@ref SiC), [`MoO₃`](@ref MoO₃), ...), we start by defining the
permittivities from scratch. 

For this we use the [`@permittivity`](@ref @permittivity) macro. This macro
expects the name of the material we want to define and a function that takes the
wavelength `λ` in meters and returns the full permittivity tensor.

A simple example would be to define Glass with a constant permittivity of 1.5:
```@example tutorial1
using LinearAlgebra

@permittivity "Glass" λ -> Diagonal(ones(3)) * 1.5
```
We first load the `LinearAlgebra` standard library to have access to the
function `Diagonal`, so we don't have to write out the ``3 \times 3`` diagonal
matrix ourselves. In this example we passed a anonymous function (`->` syntax)
to the permittivity macro, but we could have also defined the function
beforehand:
```julia
f(λ) = Diagonal(ones(3)) * 1.5

@permittivity "Glass" f
```

The macro has defined a function called "Glass," which returns a structure of
type [`Layer`](@ref Layer). Additionally, it has created a function `ϵ_Glass`,
which is identical to the input function we provided:

```@example tutorial1
Glass()
```

To calculate the reflection of a Glass--Air interface, we must define an Air
layer. To accomplish this, we can simply invoke the [`Layer`](@ref Layer)
function without providing any arguments. This will create a Layer with a
constant permittivity of 1:
```@example tutorial1
Air = Layer()
```

With these two layers, we can now build up a [`LayeredStructure`](@ref
LayeredStructure):
```@example tutorial1
Interface = LayeredStructure(
    superstrate = Glass(),
    substrate = Air
)
```

!!! tip "Intermediate layers"
    More intermediate layers can be specified using the `layers` keyword, e.g.:
    ```@example tutorial1
    LayeredStructure(
        superstrate = Air,
        layers = [Glass()],
        substrate = Air
    )
    ```

!!! danger "Isotropic superstrate and substrate"
    For this code to work the superstrate and substrate need to be
    isotropic. Most of the time, adding an extra air layer solves related
    problems. 

#### Calculating optical properties

For a given wavelength `λ` and angle of incidence `α`:
```@example tutorial1
λ = 1.55e-6 # [m]
α = deg2rad(10) # [rad]
nothing #hide
```
we can calculate the [reflection](@ref calculate_reflection):
```@example tutorial1
Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ = calculate_reflection(λ, α, Interface)
```
and [transmission](@ref calculate_transmission):
```@example tutorial1
Tₚₚ, Tₛₛ, Tₚₛ, Tₛₚ = calculate_transmission(λ, α, Interface)
```
which are returned as tuples of the form `(Xₚₚ, Xₛₛ, Xₚₛ, Xₛₚ)`, where the right
index denotes the polarzation of the incident light and the left index the
polarization of the outgoing light. 

!!! tip "Circular basis"
    The output can be changed to a circular basis using the `basis` keyword
    argument: 
    ```@example tutorial1
    calculate_transmission(λ, α, Interface; basis=:circular)
    ```
    Here, the returned tuple has the form `(T_RR, T_LL, T_RL, T_LR)`.

#### Angular dependence

We can now calculate the angular dependence of the reflection of our interface
by replacing the incident angle `α` by a list of angles:
```@example tutorial1
λ = 1.55e-6 # [m]
α = deg2rad.(0:0.1:89) # [rad]

R = calculate_reflection.(λ, α, Interface)
nothing #hide
```

We now have a list of tupels `R`, where each element has the shape `(Rₚₚ, Rₛₛ,
Rₚₛ, Rₛₚ)`. If we want to split them into separate lists, we can use the
[Unzip.jl](https://github.com/bramtayl/Unzip.jl) package:
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
savefig("interface-plot.svg"); 
nothing # hide
```

![](interface-plot.svg)
