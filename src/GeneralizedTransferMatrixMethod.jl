module GeneralizedTransferMatrixMethod

# Types.jl
export
    Layer,
    LayeredStructure

# OpticsFunctions.jl
export
    reflection_coeffs,
    transmission_coeffs,
    calculate_reflection,
    calculate_transmission

# Permitivities.jl
export
    @permittivity,
    ϵ_vacuum,
    Ag, ϵ_Ag,
    Au, ϵ_Au,
    SiC, ϵ_SiC,
    MoO₃, ϵ_MoO₃

using LinearAlgebra
using StaticArrays
using UnPack


# Physical constants
"Speed of light in vacuum ``[\\frac{m}{s}]``."
const c₀ = 299792458
"Vacuum permittivity ``[\\frac{F}{m}]``."
const ϵ₀ = 8.8541878188e-12
"Vacuum permeability ``[\\frac{N}{A^2}]``."
const μ₀ = 1.25663706127e-6

# References
"Dictionary of references used in this Package."
const References = Dict(
    "MacKay" => "Mackay, T. G. & Lakhtakia, A. The Transfer-Matrix Method in
Electromagnetics and Optics. vol. 1 (2020).",
    "Byrnes" => "Byrnes, S. J. Multilayer Optical
Calculations. arXiv:1603.02720 [Physics], December 30,
2020. http://arxiv.org/abs/1603.02720.",
    "Yeh" => "Yeh, P. Optical Waves in Layered Media. Wiley Series in Pure
and Applied Optics. Wiley, 2005.",
    "Álvarez-Pérez" => "Álvarez-Pérez, G. et al. Infrared Permittivity of the
Biaxial van Der Waals Semiconductor α-MoO₃ from Near- and Far-Field Correlative
Studies. Adv. Mater. 32, 1908176  (2020)."
)

 # Transfer matrix core
include("TMM.jl")
# Calculate optical quantities from transfer matrix
include("OpticsFunctions.jl")
# Composite types (Layer, LayeredStructure)
include("Types.jl")
# @permittivity macro and some predefined permittivities
include("Permittivities.jl")
# Helpers
include("HelperFunctions.jl")

end # module
