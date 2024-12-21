module GeneralizedTransferMatrixMethod

# Types
export
    Layer,
    LayeredStructure,
    LayerProperties,
    StructureProperties

# Functions
export
    calculate_layer_properties,
    calculate_structure_properties,
    reflection,
    transmission,
    reflection_coeffs,
    transmission_coeffs

# Permitivities
export
    @permittivity,
    ϵ_vacuum,
    Ag, ϵ_Ag,
    Au, ϵ_Au,
    SiC, ϵ_SiC,
    MoO₃, ϵ_MoO₃

using Reexport
using LinearAlgebra
using DelimitedFiles
using StaticArrays

@reexport using Unitful


# Physical constants
"Speed of light in vacuum in ``\\frac{m}{s}``."
const c₀ = 299792458
"Vacuum permittivity in ``\\frac{F}{m}``."
const ϵ₀ = 8.8541878188e-12
"Vacuum permeability in ``\\frac{N}{A^2}``."
const μ₀ = 1.25663706127e-6


include("Permittivities.jl")

include("Types.jl")
include("Matrices.jl")
include("Functions.jl")

end # module
