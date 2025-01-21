##------------------------------------------------------------------------------
## Permittivity macro
##
## This macro is used to define all permittivities!
##------------------------------------------------------------------------------

"""
    permittivity(name, func)

This macro takes the `name` of a material and the permittivity `func`,
depending on λ and returning the permittivity matrix. It generates the
permittivity funcion `ϵ_Name` and the Layer function `Name`.

# Note

The first letter of the name is always capitalised.
"""
macro permittivity(name, func)
    local name_u = uppercasefirst(name)

    local name_sym = Symbol(name_u)
    local eps = Symbol("ϵ_" * name_u)
    quote
        # Define permittivity function
        """
            $($(esc(eps)))(λ)

        Calculate permittivity tensor of $($(esc(name_u))).

        # Arguments
        - `λ`: Wavelength ``[m]``.
        """
        function $(esc(eps))(λ)
            λ = convert_to_wavelength(λ) # HelperFunctions.jl
            f = parse_permittivity($(esc(func)))
            f(λ)
        end

        # Define Layer with keyword arguments
        """
            $($(esc(name_sym)))(;d = 0, θ = 0, ϕ = 0, ψ = 0)

        Define a `Layer` of $($(esc(name_u))) using keyword arguments.

        # Arguments
        Suitable quantities with units from `Unitful` also work.
        - `d`: Thickness of the layer ``[m]`` (default: 0).
        - `θ`: θ Euler angle ``[rad]`` (default: 0).
        - `ϕ`: ϕ Euler angle ``[rad]`` (default: 0).
        - `ψ`: ψ Euler angle ``[rad]`` (default: 0).
        """
        function $(esc(name_sym))(;d = 0, θ = 0, ϕ = 0, ψ = 0)
            $(esc(name_sym))(d, θ, ϕ, ψ)
        end
        # General Function
        """
            $($(esc(name_sym)))(d, θ, ϕ, ψ)

        Define a `Layer` of $($(esc(name_u))).

        # Arguments
        - `d`: Thickness of the layer ``[m]``.
        - `θ`: θ Euler angle ``[rad]``.
        - `ϕ`: ϕ Euler angle ``[rad]``.
        - `ψ`: ψ Euler angle ``[rad]``.
        """
        function $(esc(name_sym))(d, θ, ϕ, ψ)
            Layer(
                ϵ = $(esc(eps)),
                d = d,
                θ = θ,
                ϕ = ϕ,
                ψ = ψ
            )
        end
    end
end
"""
    parse_permittivity(func)

Function used inside the `@permittivity` macro to dispatch on for
extensions. The default implementation just returns the original function.
"""
parse_permittivity(func) = func

##------------------------------------------------------------------------------
## Helper functions
##------------------------------------------------------------------------------

"""
    lorentz_osc(f, fₗₒ, fₜₒ, γ)

Calculate a single lorentz oscillator (without ϵ∞). The definition is taken from [^1].

# Arguments

- `f`: Frequency.
- `fₗₒ`: Frequency of the longitudinal optical phonon.
- `fₜₒ`: Frequency of the transverse optical phonon.
- `γ`: Damping factor of the Lorentzian line shape.

# References

[^1]: $(References["Álvarez-Pérez"])
"""
function lorentz_osc(f, fₗₒ, fₜₒ, γ)
    numerator = fₗₒ^2 - f^2 - 1im * γ * f
    denominator = fₜₒ^2 - f^2 - 1im * γ * f
    return numerator / denominator
end


"""
    ϵ_drude(f, fₚ, γ)

Calculate the permitivity from Drude model.

# Arguments

- `f`: Frequency.
- `fₚ`: Plasma frequency.
- `γ`: Mean collision rate.
"""
function ϵ_drude(f, fₚ, γ, ϵ∞ = 1.0)
   ϵ∞ - fₚ^2 / (f^2 + 1im * f * γ)
end
##------------------------------------------------------------------------------
## Vacuum
##------------------------------------------------------------------------------

"""
    ϵ_vacuum(λ)

Calculate the relative permittivity of vacuum. Always returns the identity
matrix.
"""
ϵ_vacuum(λ) = SIdentity

"""
    μ_vacuum(λ)

Calculate the relative permeability of vacuum. Always returns the identity
matrix.
"""
μ_vacuum(λ) = SIdentity

"Optical rotation tensor"
ξ_vacuum(λ) = @SMatrix zeros(3,3)
"Optical rotation tensor"
ζ_vacuum(λ) = @SMatrix zeros(3,3)

##------------------------------------------------------------------------------
## MODELLED PERMITIVITIES
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
## Metals
##------------------------------------------------------------------------------

@permittivity "Ag" λ -> let f = c₀ / λ,
    fₚ_Ag = 2.152e15,
    γ_Ag = 1/17e-15,
    ϵ∞_Ag = 5;

    SIdentity .* ϵ_drude(f, fₚ_Ag, γ_Ag, ϵ∞_Ag)
end

@permittivity "Au" λ -> let f = c₀ / λ,
    fₚ_Au = 2.183e15,
    γ_Au = 1.7410e13,
    ϵ∞_Au = 9.84;

    SIdentity .* ϵ_drude(f, fₚ_Au, γ_Au, ϵ∞_Au)
end

##------------------------------------------------------------------------------
## Polar materials
##------------------------------------------------------------------------------

# Reference: https://doi.org/10.1364/OL.34.002667
@permittivity "SiC" λ -> let f = 1 / (λ * 1e2),
    # Parameters from paper
    ϵ∞ = 6.5,
    fₗₒ = 972, # [cm⁻¹]
    fₜₒ = 796, # [cm⁻¹]
    γ = 3.75; # [cm⁻¹]

    SIdentity .* ϵ∞ * lorentz_osc(f, fₗₒ, fₜₒ, γ)
end

##------------------------------------------------------------------------------
## Hyperbolic materials
##------------------------------------------------------------------------------


# Reference: https://doi.org/10.1002/adma.201908176
@permittivity "MoO₃" λ -> Diagonal(@SVector [ϵ_x_MoO₃(λ), ϵ_y_MoO₃(λ), ϵ_z_MoO₃(λ)])

"""
    ϵ_x_MoO₃(λ)

Calculate the x principal component of the permitivity tensor of MoO₃. The
parameters are taken from [^1].

# Arguments

- `λ`: Wavelength `[m]`.

# References

[^1]: $(References["Álvarez-Pérez"])
"""
function ϵ_x_MoO₃(λ)
    # Convert λ in meter to frequency in cm⁻¹
    f = 1 / (λ * 1e2)
    # Parameters from paper
    ϵ∞_x = 5.78
    fₗₒ_x₁, fₗₒ_x₂, fₗₒ_x₃ = 534.3, 963.0, 999.2 # [cm⁻¹]
    fₜₒ_x₁, fₜₒ_x₂, fₜₒ_x₃ = 506.7, 821.4, 998.7 # [cm⁻¹]
    γ_x₁, γ_x₂, γ_x₃ = 49.1, 6.0, 0.35# [cm⁻¹]

    return (ϵ∞_x * lorentz_osc(f, fₗₒ_x₁, fₜₒ_x₁, γ_x₁)
            * lorentz_osc(f, fₗₒ_x₂, fₜₒ_x₂, γ_x₂)
            * lorentz_osc(f, fₗₒ_x₃, fₜₒ_x₃, γ_x₃))
end


"""
    ϵ_y_MoO₃(λ)

Calculate the y principal component of the permitivity tensor of MoO₃. The
parameters are taken from [^1].

# Arguments

- `λ`: Wavelength `[m]`.

# References

[^1]: $(References["Álvarez-Pérez"])
"""
function ϵ_y_MoO₃(λ)
    # Convert λ in meter to frequency in cm⁻¹
    f = 1 / (λ * 1e2)
    # Parameters from paper
    ϵ∞_y = 6.07
    fₗₒ_y = 850.1 # [cm⁻¹]
    fₜₒ_y = 544.6 # [cm⁻¹]
    γ_y = 9.5 # [cm⁻¹]

    return ϵ∞_y * lorentz_osc(f, fₗₒ_y, fₜₒ_y, γ_y)
end


"""
    ϵ_z_MoO₃(λ)

Calculate the z principal component of the permitivity tensor of MoO₃. The
parameters are taken from [^1].

# Arguments

- `λ`: Wavelength `[m]`.

# References

[^1]: $(References["Álvarez-Pérez"])
"""
function ϵ_z_MoO₃(λ)
    # Convert λ in meter to frequency in cm⁻¹
    f = 1 / (λ * 1e2)
    # Parameters from paper
    ϵ∞_z = 4.47
    fₗₒ_z = 1006.9 # [cm⁻¹]
    fₜₒ_z = 956.7 # [cm⁻¹]
    γ_z = 1.5 # [cm⁻¹]

    return ϵ∞_z * lorentz_osc(f, fₗₒ_z, fₜₒ_z, γ_z)
end
