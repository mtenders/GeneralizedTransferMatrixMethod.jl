module RefractiveIndexExt

using LinearAlgebra
using StaticArrays

import GeneralizedTransferMatrixMethod as TMM

using RefractiveIndex


# Helpers
"""
    to_permittivity(m::RefractiveMaterial, λ)

Calculate permittivity from `RefractiveMaterial` with `λ` in ``[m]``.
"""
function to_permittivity(m::RefractiveMaterial, λ)
    λ_μm = λ * 1e6
    m.dispersion(λ_μm)^2
end

function to_permittivity(m::RefractiveMaterial{RefractiveIndex.TabulatedNK}, λ)
    λ_μm = λ * 1e6
    (m.dispersion.n(λ_μm) + im * m.dispersion.k(λ_μm))^2
end

function to_permittivity(m::RefractiveMaterial{RefractiveIndex.TabulatedK}, λ)
    throw(ArgumentError("`RefractiveMaterial{TabulatedK}` is not supported."))
end

# Dispatch
function TMM.parse_permittivity(m::RefractiveMaterial)
    λ -> to_permittivity(m, λ) * TMM.SIdentity
end
function TMM.parse_permittivity(m::Vector{T}) where T<:RefractiveMaterial
    if length(m) == 1
        TMM.parse_permittivity(only(m))
    elseif length(m) == 2
        λ -> Diagonal(@SVector [to_permittivity(m[1], λ), to_permittivity(m[2], λ),
                                to_permittivity(m[2], λ)])
    elseif length(m) == 3
        λ -> Diagonal(@SVector [to_permittivity(m[1], λ), to_permittivity(m[2], λ),
                                to_permittivity(m[3], λ)])
    else
        throw(ArgumentError("The array of refractive materials has the wrong dimension!"))
    end
end

end # module
