module UnitfulExt

import GeneralizedTransferMatrixMethod as TMM

using Unitful

# Convert common units to wavelength
convert_to_wavelength(x::Unitful.Length) = ustrip(u"m", x)
convert_to_wavelength(x::Quantity{<:Any,Unitful.𝐋^-1}) = ustrip(u"m", 1/x)
convert_to_wavelength(x::Quantity{<:Any,Unitful.𝐓^-1}) = ustrip(u"m", 1/x *
    Unitful.c0)
convert_to_wavelength(x::Quantity{<:Any,Unitful.𝐋^2 * Unitful.𝐌 * Unitful.𝐓^-2}) =
    ustrip(u"m", Unitful.h / x * Unitful.c0)


function TMM.calculate_reflection(λ::Quantity, α, strct; basis=:linear)
    λ_SI = convert_to_wavelength(λ)
    TMM.calculate_reflection(λ_SI, α, strct; basis=basis)
end

function TMM.calculate_transmission(λ::Quantity, α, strct; basis=:linear)
    λ_SI = convert_to_wavelength(λ)
    TMM.calculate_transmission(λ_SI, α, strct; basis=basis)
end

# Dispatch on Unitful.Length (angles work out of the box)
function TMM.Layer(ϵ::Function, μ::Function, ξ::Function, ζ::Function,
                   d::Unitful.Length, θ, ϕ, ψ)
    TMM.Layer(ϵ, μ, ξ, ζ, ustrip(u"m",d), θ, ϕ, ψ)
end

end
