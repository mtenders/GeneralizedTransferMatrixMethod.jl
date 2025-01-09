module UnitfulExt

using Unitful: ustrip

using GeneralizedTransferMatrixMethod

"""
    calculate_reflection(λ::Unitful.Length, α, strct; basis=:linear)

Calculate the reflection of a layered structure.

# Arguments

- `λ::Unitful.Length`: Wavelength.
- `α`: Angle of incidence. (`u"°"` doesn't have any dimension)
- `strct`: Structure of type `LayeredStructure`.
- `basis=:linear`: Returns `(Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ)` in `:linear` basis and `(R_RR,
  R_LL, R_RL, R_LR)` in `:circular` basis.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), layers=[MoO₃(d = 1e-6)],
substrate = Layer());
julia> calculate_reflection(12.5u"μm", 23u"°", S)
(0.19265761641397677, 0.31745983790535026, 0.0, 0.0)
```
"""
function calculate_reflection(λ::Unitful.Length, α, strct; basis=:linear)
    calculate_reflection(ustrip(u"m",λ), α, strct, basis=basis)
end

"""
    calculate_transmission(λ::Unitful.Length, α, strct; basis=:linear)

Calculate the transmission of a layered structure.

# Arguments

- `λ::Unitful.Length`: Wavelength.
- `α`: Angle of incidence. (`u"°"` doesn't have any dimension)
- `strct`: Structure of type `LayeredStructure`.
- `basis=:linear`: Returns `(Tₚₚ, Tₛₛ, Tₚₛ, Tₛₚ)` in `:linear` basis and `(T_RR,
  T_LL, T_RL, T_LR)` in `:circular` basis.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), layers=[MoO₃(d = 1e-6)],
substrate = Layer());
julia> calculate_transmission(12.5u"μm", 23u"°", S)
(0.31887703646259946, 0.6139558802434743, 0.0, 0.0)
```
"""
function calculate_transmission(λ::Unitful.Length, α, strct; basis=:linear)
    calculate_transmission(ustrip(u"m",λ), α, strct, basis=basis)
end

end
