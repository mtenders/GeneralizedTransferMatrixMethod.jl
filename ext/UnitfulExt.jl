module UnitfulExt

import GeneralizedTransferMatrixMethod as TMM

using Unitful

# Convert common units to wavelength
TMM.convert_to_wavelength(x::Unitful.Length) = ustrip(u"m", x)
TMM.convert_to_wavelength(x::Quantity{<:Any,Unitful.𝐋^-1}) = ustrip(u"m", 1/x)
TMM.convert_to_wavelength(x::Quantity{<:Any,Unitful.𝐓^-1}) = ustrip(u"m", 1/x *
    Unitful.c0)
TMM.convert_to_wavelength(x::Quantity{<:Any,Unitful.𝐋^2 * Unitful.𝐌 *
    Unitful.𝐓^-2}) = ustrip(u"m", Unitful.h / x * Unitful.c0)

# Dispatch on Unitful.Length (angles work out of the box)
function TMM.Layer(ϵ::Function, μ::Function, ξ::Function, ζ::Function,
                   d::Unitful.Length, θ, ϕ, ψ)
    TMM.Layer(ϵ, μ, ξ, ζ, ustrip(u"m",d), θ, ϕ, ψ)
end

end
