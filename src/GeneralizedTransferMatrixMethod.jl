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
using Parameters

@reexport using Unitful

const c₀ = 299792458

include("Permittivities.jl")

include("Types.jl")
include("Matrices.jl")
include("Functions.jl")

end # module
