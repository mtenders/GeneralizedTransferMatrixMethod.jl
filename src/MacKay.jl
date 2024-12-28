module MacKay

using LinearAlgebra
using StaticArrays
using UnPack
using Unitful

using ..GeneralizedTransferMatrixMethod: ϵ₀, μ₀, c₀, euler_mat, basis_change,
    Layer

"Dictionary of references used in this Package."
const References = Dict(
    "MacKay" => "Mackay, T. G. & Lakhtakia, A. The Transfer-Matrix Method in
    Electromagnetics and Optics. vol. 1 (2020)."
)

"""
All-ones matrix as defined in Eq. 2.50 in [^1].

# References
[^1]: $(References["MacKay"])
"""
const J = @SMatrix ones(4,4)

@doc raw"""
    ν_mats(qₓ, ϵ, μ, ξ, ζ)

Calculate the coefficients of the auxiliary phasors ``e_z`` and ``h_z``, as
described in Eq. 2.46 in [^1]. Returns two matrices of the form:

```math
\begin{pmatrix}
\nu^{ee}_{zx} & & & \\
& \nu^{ee}_{zy} & & \\
& & \nu^{eh}_{zx} & \\
& & & \nu^{eh}_{zy}
\end{pmatrix}
```
and
```math
\begin{pmatrix}
\nu^{he}_{zx} & & & \\
& \nu^{he}_{zy} & & \\
& & \nu^{hh}_{zx} & \\
& & & \nu^{hh}_{zy}
\end{pmatrix}.
```

# Arguments

- `qₓ`: Normalized in-plane component of the wavevector.
- `ϵ`: Permittivity tensor.
- `μ`: Permeability tensor.
- `ξ`, `ζ`: Optical rotation tensors.

# Examples
```jldoctest
julia> ν_mats(0.3, 3ones(3,3), ones(3,3), zeros(3,3), zeros(3,3))
([-1.0 0.0 0.0 0.0;
0.0 -1.0 0.0 0.0;
0.0 0.0 0.0 0.0;
0.0 0.0 0.0 -3.3356409519815207e-10],
[0.0 0.0 0.0 0.0;
0.0 1.0006922855944562e-9 0.0 0.0;
0.0 0.0 -1.0 0.0;
0.0 0.0 0.0 -1.0])
```

# References
[^1]:""" * References["MacKay"]
function ν_mats(qₓ, ϵ, μ, ξ, ζ)
    denom = ϵ[3,3] * μ[3,3] - ξ[3,3] * ζ[3,3]

    ν_ee_zx = - (μ[3,3] * ϵ[3,1] - ξ[3,3] *  ζ[3,1]           ) / denom
    ν_ee_zy = - (μ[3,3] * ϵ[3,2] - ξ[3,3] * (ζ[3,2] - qₓ / c₀)) / denom
    ν_eh_zx =   (ξ[3,3] * μ[3,1] - μ[3,3] *  ξ[3,1]           ) / denom
    ν_eh_zy =   (ξ[3,3] * μ[3,2] - μ[3,3] * (ξ[3,2] + qₓ / c₀)) / denom
    ν_he_zx =   (ζ[3,3] * ϵ[3,1] - ϵ[3,3] *  ζ[3,1]           ) / denom
    ν_he_zy =   (ζ[3,3] * ϵ[3,1] - ϵ[3,3] * (ζ[3,2] - qₓ / c₀)) / denom
    ν_hh_zx = - (ϵ[3,3] * μ[3,1] - ζ[3,3] *  ξ[3,1]           ) / denom
    ν_hh_zy = - (ϵ[3,3] * μ[3,2] - ζ[3,3] * (ξ[3,2] + qₓ / c₀)) / denom

    ν_e = Diagonal(@SVector [ν_ee_zx, ν_ee_zy, ν_eh_zx, ν_eh_zy])
    ν_h = Diagonal(@SVector [ν_he_zx, ν_he_zy, ν_hh_zx, ν_hh_zy])

    ν_e, ν_h
end

"""
    P_mat(ω, qₓ, ϵ, μ, ξ, ζ)

Calculate the ``P`` matrix according to Eq. 2.49 in [^1].

# Arguments

- `ω`: Angular frequency `[1/s]` ``[\\frac{1}{s}]``.
- `qₓ`: Normalized in-plane component of the wavevector.
- `ϵ`: Permittivity tensor.
- `μ`: Permeability tensor.
- `ξ`, `ζ`: Optical rotation tensors.

# Examples
```jldoctest
julia> P_mat(1.8e14, 0.3, 3ones(3,3), ones(3,3), zeros(3,3), zeros(3,3))
4×4 SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 -1.80125e5   0.0           0.0         0.0
  0.0        -1.80125e5     0.0         0.0
  0.0         0.000180249  -1.80125e5   0.0
  0.0         0.0           0.0        -1.80125e5
```

# References
[^1]: $(References["MacKay"])
"""
function P_mat(ω, qₓ, ϵ, μ, ξ, ζ)
    ν_e, ν_h = ν_mats(qₓ, ϵ, μ, ξ, ζ)

    # Matrices in the 3 terms of the sum
    P₁ = @SMatrix [
	ζ[2, 1]   ζ[2, 2]  μ[2, 1]  μ[2, 2]
        -ζ[1, 1]  -ζ[1, 2] -μ[1, 1] -μ[1, 2]
	-ϵ[2, 1]  -ϵ[2, 2] -ξ[2, 1] -ξ[2, 2]
	ϵ[1, 1]   ϵ[1, 2]  ξ[1, 1]  ξ[1, 2]
    ]
    P₂ = Diagonal(@SVector [ζ[2,3] + qₓ / c₀, -ζ[1,3], -ϵ[2,3], ϵ[1,3]])
    P₃ = Diagonal(@SVector [μ[2,3], -μ[1,3], ξ[2,3] + qₓ / c₀, ξ[1,3]])

    # P
    ω * (P₁ + P₂ * J * ν_e + P₃ * J * ν_h)
end

"""
    M_mat(d, ω, qₓ, ϵ, μ, ξ, ζ)

Calculate the ``M`` matrix according to Eq. 2.52 in [^1].

# Arguments

- `d`: Thickness ``[m]``.
- `ω`: Angular frequency ``[\\frac{1}{s}]``.
- `qₓ`: Normalized in-plane component of the wavevector.
- `ϵ`: Permittivity tensor.
- `μ`: Permeability tensor.
- `ξ`, `ζ`: Optical rotation tensors.

# Examples
```jldoctest
julia> M_mat(1e-6, 1.8e14, 0.3, 3ones(3,3), ones(3,3), zeros(3,3), zeros(3,3))
4×4 SMatrix{4, 4, ComplexF64, 16} with indices SOneTo(4)×SOneTo(4):
 0.983821+0.179152im          0.0-0.0im               0.0-0.0im            0.0-0.0im
      0.0-0.0im          0.983821+0.179152im          0.0-0.0im            0.0-0.0im
      0.0-0.0im       3.22921e-11-1.77333e-10im  0.983821+0.179152im       0.0-0.0im
      0.0-0.0im               0.0-0.0im               0.0-0.0im       0.983821+0.179152im
```

# References
[^1]: $(References["MacKay"])
"""
function M_mat(d, ω, qₓ, ϵ, μ, ξ, ζ)
    P = P_mat(ω, qₓ, ϵ, μ, ξ, ζ)

    exp(-im * P * d)
end

"""
    K_mat(n, θ)

Calculate the ``K`` matrix according to Eq. 3.16 in [^1].

# Arguments

- `n`: Refractive index.
- `θ`: Angle between the k-vector and the surface normal.

# Examples
```jldoctest
julia> K_mat(1.8, π/4)
4×4 SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 0.0          0.707107     0.0         -0.707107
 1.0          0.0          1.0          0.0
 0.00337852   0.0         -0.00337852   0.0
 0.0         -0.00477795   0.0         -0.00477795
```

# References
[^1]: $(References["MacKay"])
"""
function K_mat(n, θ)
    N = n * sqrt(ϵ₀ / μ₀)

    @SMatrix [
	0           cos(θ)  0          -cos(θ)
	1           0       1           0
	N * cos(θ)  0      -N * cos(θ)  0
	0          -N       0          -N
    ]
end

"""
    calculate_layer_M(λ, qₓ, layer)

Calculate the ``M`` matrix for a single layer.

# Arguments
- `λ`: Wavelength ``[m]``.
- `qₓ`: Normalized in-plane component of the wavevector.
- `layer`: Material `Layer`.

# Examples
```jldoctest
julia> calculate_layer_M(12.3e-6, 0.5, Layer(d = 2e-6))
4×4 StaticArraysCore.SMatrix{4, 4, ComplexF64, 16} with indices SOneTo(4)×SOneTo(4):
 0.63346+0.0im             0.0+0.0im             0.0+0.0im          0.0-252.451im
     0.0+0.0im         0.63346+0.0im             0.0+336.601im      0.0+0.0im
     0.0+0.0im             0.0+0.00177875im  0.63346+0.0im          0.0+0.0im
     0.0-0.00237167im      0.0+0.0im             0.0+0.0im      0.63346+0.0im
```
"""
function calculate_layer_M(λ, qₓ, layer)
    ω = 2π * c₀ / λ
    eul = euler_mat(layer)

    ϵ = ϵ₀ * basis_change(layer.ϵ(λ), inv(eul)) # TODO double check the
    # inversion is correct, something is still weird with theta and psi
    μ = μ₀ * basis_change(layer.μ(λ), inv(eul))
    ξ = basis_change(layer.ξ(λ), inv(eul))
    ζ = basis_change(layer.ζ(λ), inv(eul))

    M_mat(layer.d, ω, qₓ, ϵ, μ, ξ, ζ)
end

"""
    calculate_α_out(λ, α, strct)

Calculate the outgoing angle for a layered structure `strct` and angle of
incidence `α`. Returns the tuple `(α_out, n_inc, n_out, qₓ)`, where `α_out` is the
outgoing angle, `n_inc` is the refractive index of the incident medium, `n_out` is
the refractive index of the outgoing medium, and `qₓ` is the normalized in-plane
wavevector component.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), substrate = SiC());
julia> calculate_α_out(12e-6, deg2rad(14), S)
(0.0014914047759539242 - 0.04676247666123382im, 1.0 + 0.0im, 0.1648892944698832
+ 5.166277293025589im, 0.24192189559966773 + 0.0im)
```
"""
function calculate_α_out(λ, α, strct)
    # Only works for isotropic in- and outgoing media
    @unpack superstrate, layers, substrate = strct

    n_inc = sqrt(Complex(superstrate.ϵ(λ)[1,1] * superstrate.μ(λ)[1,1]))
    n_out = sqrt(Complex(substrate.ϵ(λ)[1,1] * substrate.μ(λ)[1,1]))
    qₓ = n_inc * sin(α)

    α_out = asin(Complex(qₓ / n_out))

    α_out, n_inc, n_out, qₓ
end


"""
    calculate_structure_M(λ, α, strct)

Calculate the ``M`` matrix for a layered structure.

# Arguments
- `λ`: Wavelength ``[m]``.
- `α`: Angle of incidence.
- `strct`: Structure of type `LayeredStructure`.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), substrate=Au());
julia> calculate_structure_M(12.3e-6, deg2rad(20), S)
4×4 StaticArraysCore.SMatrix{4, 4, ComplexF64, 16} with indices SOneTo(4)×SOneTo(4):
 0.501778-0.00554054im       0.0+0.0im         0.498222+0.00554054im       0.0+0.0im
      0.0+0.0im         0.471735-0.00589861im       0.0+0.0im         -0.46795-0.0058937im
 0.498222+0.00554054im       0.0+0.0im         0.501778-0.00554054im       0.0+0.0im
      0.0+0.0im         -0.46795-0.0058937im        0.0+0.0im         0.471735-0.00589861im
```
"""
function calculate_structure_M(λ, α, strct)
    @unpack superstrate, layers, substrate = strct

    α_out, n_inc, n_out, qₓ = calculate_α_out(λ, α, strct)

    K_sup = K_mat(n_inc, α)
    K_sub = K_mat(n_out, α_out)

    if !isempty(layers)
	lay_Ms = calculate_layer_M.(λ, qₓ, layers)

	M = inv(K_sub) * prod(reverse(lay_Ms)) * K_sup
    else
	M = inv(K_sub) * K_sup
    end

    M
end

# Basis change stuff
"""
Basis tarnsformation matrix to change from a linear polarization basis to a
circular polarization basis.
"""
const T_mat = 1/sqrt(2) * [
	1 -im
	1 im
]

"""
    circ_coeffs(cₚₚ, cₛₛ, cₚₛ, cₛₚ)

Transform transmission/reflection coefficients from a linear to a circula
basis. The output has the form `(c_RR, c_LL, c_RL, c_LR)`.

# Examples
```jldoctest
julia> circ_coeffs(0.23, 0.24, 0.14, 0.14)
(0.235 + 0.0im, 0.235 + 0.0im, -0.0050000000000000044 + 0.14im,
-0.0050000000000000044 - 0.14im)
```
"""
function circ_coeffs(cₚₚ, cₛₛ, cₚₛ, cₛₚ)
    c_linear = [
        cₚₚ cₚₛ
        cₛₚ cₛₛ
    ]

    c_RR, c_RL, c_LR, c_LL = basis_change(c_linear, inv(T_mat))
    return c_RR, c_LL, c_RL, c_LR
end

"""
    basis_selector(coeffs, basis)

Return linear input `coeffs` in given `basis`, either `:linear` or `:circular`.
"""
function basis_selector(coeffs, basis)
    if basis == :linear
        coeffs
    elseif basis == :circular
        circ_coeffs(coeffs...)
    else
        throw(ArgumentError("$(repr(basis)) basis not defined. Only :linear and :circular
bases are defined."))
    end
end

# Calculate coefficients
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
(0.9916655707450739 + 0.02489351332010354im, -0.9926693143906693 -
0.02200268526284087im, 0.0 + 0.0im, 0.0 + 0.0im)
```
"""
function reflection_coeffs(λ, α, strct; basis=:linear)
    M = calculate_structure_M(λ, α, strct)

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
julia> M = calculate_structure_M(12.3e-6, deg2rad(20), S);
julia> reflection_coeffs(M)
(0.9916655707450739 + 0.02489351332010354im, -0.9926693143906693 -
0.02200268526284087im, 0.0 + 0.0im, 0.0 + 0.0im)
```
"""
function reflection_coeffs(M; basis=:linear)
    denom = M[4,4] * M[3,3] - M[4,3] * M[3,4]

    rₚₚ = (M[3,2] * M[4,3] - M[3,3] * M[4,2]) / denom
    rₛₛ = (M[4,1] * M[3,4] - M[3,1] * M[4,4]) / denom
    rₚₛ = (M[4,2] * M[3,4] - M[4,4] * M[3,2]) / denom
    rₛₚ = (M[3,1] * M[4,3] - M[3,3] * M[4,1]) / denom

    basis_selector((rₚₚ, rₛₛ, rₚₛ, rₛₚ), basis)
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
(0.08647336439281672 - 0.36672675709549724im, 0.07570669303896788 -
0.35120281878473im, 0.0 + 0.0im, 0.0 + 0.0im)
```
"""
function transmission_coeffs(λ, α, strct; basis=:linear)
    M = calculate_structure_M(λ, α, strct)

    transmission_coeffs(M; basis=basis)
end
"""
    transmission_coeffs(M; basis=:linear)

Calculate transmission coefficients from `M` matrix. Returns `(tₚₚ, tₛₛ, tₚₛ,
tₛₚ)` in `:linear` basis and `(t_RR, t_LL, t_RL, t_LR)` in `:circular` basis.

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), substrate = SiC());
julia> M = calculate_structure_M(12e-6, deg2rad(18), S);
julia> transmission_coeffs(M)
(0.08647336439281672 - 0.36672675709549724im, 0.07570669303896788 -
0.35120281878473im, 0.0 + 0.0im, 0.0 + 0.0im)
```
"""
function transmission_coeffs(M; basis=:linear)
    denom = M[4,4] * M[3,3] - M[4,3] * M[3,4]
    rₚₚ, rₛₛ, rₚₛ, rₛₚ = reflection_coeffs(M)

    tₛₛ = M[1,1] + M[1,3] * rₛₛ + M[1,4] * rₛₚ
    tₛₚ = M[2,1] + M[2,3] * rₛₛ + M[2,4] * rₛₚ
    tₚₛ = M[1,2] + M[1,3] * rₚₛ + M[1,4] * rₚₚ
    tₚₚ = M[2,2] + M[2,3] * rₚₛ + M[2,4] * rₚₚ

    basis_selector((tₚₚ, tₛₛ, tₚₛ, tₛₚ), basis)
end

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
julia> S = LayeredStructure(superstrate=Layer(), layers=[MoO₃(d = 1e-6)],
substrate = Layer());
julia> calculate_reflection(12.5e-6, deg2rad(23), S)
(0.19265761641397677, 0.31745983790535026, 0.0, 0.0)
```
"""
function calculate_reflection(λ, α, strct; basis=:linear)
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
julia> S = LayeredStructure(superstrate=Layer(), layers=[MoO₃(d = 1e-6)],
substrate = Layer());
julia> calculate_transmission(12.5e-6, deg2rad(23), S)
(0.31887703646260473, 0.6139558802434731, 0.0, 0.0)
```
"""
function calculate_transmission(λ, α, strct; basis=:linear)
    @unpack superstrate, layers, substrate = strct

    α_out, n_inc, n_out, qₓ = calculate_α_out(λ, α, strct)

    t = transmission_coeffs(λ, α, strct; basis=basis)

    real(n_out * cos(α_out)) / real(n_inc * cos(α)) .* abs2.(t)
end

### UNITFUL
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
