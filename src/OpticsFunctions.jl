##------------------------------------------------------------------------------
## E-Field quantities
##------------------------------------------------------------------------------

"""
    reflection_coeffs(λ, α, strct; basis=:linear)

Calculate the reflection coefficients.

# Arguments

- `λ`: Wavelength ``[m]``.
- `α`: Angle of incidence.
- `strct`: Structure of type `LayeredStructure`.
- `basis=:linear`: Returns `(rₚₚ, rₛₛ, rₚₛ, rₛₚ)` in `:linear` basis and `(r_RR,
  r_LL, r_RL, r_LR)` in `:circular` basis.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), substrate=Au());

julia> reflection_coeffs(12.3e-6, deg2rad(20), S)
(0.9916655707450739 + 0.02489351332010354im, -0.9926693143906693 - 0.02200268526284087im, 0.0 + 0.0im, -0.0 + 0.0im)
```
"""
function reflection_coeffs(λ, α, strct; basis=:linear)
    λ = convert_to_wavelength(λ) # HelperFunctions.jl
    M = calculate_structure_M(λ, α, strct) # TMM.jl

    reflection_coeffs(M; basis=basis)
end
"""
    reflection_coeffs(M; basis=:linear)

Calculate reflection coefficients from `M` matrix. Returns `(rₚₚ, rₛₛ, rₚₛ,
rₛₚ)` in `:linear` basis and `(r_RR, r_LL, r_RL, r_LR)` in `:circular`
basis.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), substrate=Au());

julia> M = GeneralizedTransferMatrixMethod.calculate_structure_M(12.3e-6, deg2rad(20), S);

julia> reflection_coeffs(M)
(0.9916655707450739 + 0.02489351332010354im, -0.9926693143906693 - 0.02200268526284087im, 0.0 + 0.0im, -0.0 + 0.0im)
```
"""
function reflection_coeffs(M; basis=:linear)
    denom = M[4,4] * M[3,3] - M[4,3] * M[3,4]

    rₚₚ = (M[3,2] * M[4,3] - M[3,3] * M[4,2]) / denom
    rₛₛ = (M[4,1] * M[3,4] - M[3,1] * M[4,4]) / denom
    rₚₛ = (M[4,2] * M[3,4] - M[4,4] * M[3,2]) / denom
    rₛₚ = (M[3,1] * M[4,3] - M[3,3] * M[4,1]) / denom

    basis_selector((rₚₚ, rₛₛ, rₚₛ, rₛₚ), basis) # HelperFunctions.jl
end

"""
    transmission_coeffs(λ, α, strct; basis=:linear)

Calculate the transmission coefficients.

# Arguments

- `λ`: Wavelength ``[m]``.
- `α`: Angle of incidence.
- `strct`: Structure of type `LayeredStructure`.
- `basis=:linear`: Returns `(tₚₚ, tₛₛ, tₚₛ, tₛₚ)` in `:linear` basis and `(t_RR,
  t_LL, t_RL, t_LR)` in `:circular` basis.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), substrate = SiC());

julia> transmission_coeffs(12e-6, deg2rad(18), S)
(0.08647336439281672 - 0.3667267570954972im, 0.07570669303896788 - 0.35120281878473003im, 0.0 + 0.0im, 0.0 + 0.0im)
```
"""
function transmission_coeffs(λ, α, strct; basis=:linear)
    λ = convert_to_wavelength(λ) # HelperFunctions.jl
    M = calculate_structure_M(λ, α, strct) # TMM.jl

    transmission_coeffs(M; basis=basis)
end
"""
    transmission_coeffs(M; basis=:linear)

Calculate transmission coefficients from `M` matrix. Returns `(tₚₚ, tₛₛ, tₚₛ,
tₛₚ)` in `:linear` basis and `(t_RR, t_LL, t_RL, t_LR)` in `:circular` basis.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), substrate = SiC());

julia> M = GeneralizedTransferMatrixMethod.calculate_structure_M(12e-6, deg2rad(18), S);

julia> transmission_coeffs(M)
(0.08647336439281672 - 0.3667267570954972im, 0.07570669303896788 - 0.35120281878473003im, 0.0 + 0.0im, 0.0 + 0.0im)
```
"""
function transmission_coeffs(M; basis=:linear)
    denom = M[4,4] * M[3,3] - M[4,3] * M[3,4]
    rₚₚ, rₛₛ, rₚₛ, rₛₚ = reflection_coeffs(M)

    tₛₛ = M[1,1] + M[1,3] * rₛₛ + M[1,4] * rₛₚ
    tₛₚ = M[2,1] + M[2,3] * rₛₛ + M[2,4] * rₛₚ
    tₚₛ = M[1,2] + M[1,3] * rₚₛ + M[1,4] * rₚₚ
    tₚₚ = M[2,2] + M[2,3] * rₚₛ + M[2,4] * rₚₚ

    basis_selector((tₚₚ, tₛₛ, tₚₛ, tₛₚ), basis) # HelperFunctions.jl
end

##------------------------------------------------------------------------------
## Intensity quantities
##------------------------------------------------------------------------------

"""
    calculate_reflection(λ, α, strct; basis=:linear)

Calculate the reflection of a layered structure.

# Arguments

- `λ`: Wavelength ``[m]``.
- `α`: Angle of incidence.
- `strct`: Structure of type `LayeredStructure`.
- `basis=:linear`: Returns `(Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ)` in `:linear` basis and `(R_RR,
  R_LL, R_RL, R_LR)` in `:circular` basis.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), layers=[MoO₃(d = 1e-6)], substrate = Layer());

julia> calculate_reflection(12.5e-6, deg2rad(23), S)
(0.19265761641397677, 0.31745983790535026, 0.0, 0.0)
```
"""
function calculate_reflection(λ, α, strct; basis=:linear)
    λ = convert_to_wavelength(λ) # HelperFunctions.jl
    r = reflection_coeffs(λ, α, strct; basis=basis)
    abs2.(r)
end

"""
    calculate_transmission(λ, α, strct; basis=:linear)

Calculate the transmission of a layered structure.

# Arguments

- `λ`: Wavelength ``[m]``.
- `α`: Angle of incidence.
- `strct`: Structure of type `LayeredStructure`.
- `basis=:linear`: Returns `(Tₚₚ, Tₛₛ, Tₚₛ, Tₛₚ)` in `:linear` basis and `(T_RR,
  T_LL, T_RL, T_LR)` in `:circular` basis.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), layers=[MoO₃(d = 1e-6)], substrate = Layer());

julia> calculate_transmission(12.5e-6, deg2rad(23), S)
(0.31887703646260473, 0.6139558802434731, 0.0, 0.0)
```
"""
function calculate_transmission(λ, α, strct; basis=:linear)
    λ = convert_to_wavelength(λ) # HelperFunctions.jl
    @unpack superstrate, layers, substrate = strct

    cos_αt, n_inc, n_out, qₓ = calculate_cos_αt(λ, α, strct) # HelperFunctions.jl

    t = transmission_coeffs(λ, α, strct; basis=basis)

    real(n_out * cos_αt) / real(n_inc * cos(α)) .* abs2.(t)
end
