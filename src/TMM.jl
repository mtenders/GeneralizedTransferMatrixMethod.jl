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
    K_mat(n, cos_θ)

Calculate the ``K`` matrix according to Eq. 3.16 in [^1].

# Arguments

- `n`: Refractive index.
- `cos_θ`: Cosine of angle between the k-vector and the surface normal.

# Examples
```jldoctest
julia> K_mat(1.8, cos(π/4))
4×4 SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
 0.0          0.707107     0.0         -0.707107
 1.0          0.0          1.0          0.0
 0.00337852   0.0         -0.00337852   0.0
 0.0         -0.00477795   0.0         -0.00477795
```

# References
[^1]: $(References["MacKay"])
"""
function K_mat(n, cos_θ)
    N = n * sqrt(ϵ₀ / μ₀)

    @SMatrix [
	0           cos_θ   0          -cos_θ
	1           0       1           0
	N * cos_θ   0      -N * cos_θ   0
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
    eul = euler_mat(layer) # HelperFunctions.jl

    ϵ = ϵ₀ * basis_change(layer.ϵ(λ), inv(eul)) # HelperFunctions.jl
    μ = μ₀ * basis_change(layer.μ(λ), inv(eul)) # HelperFunctions.jl
    ξ = basis_change(layer.ξ(λ), inv(eul)) # HelperFunctions.jl
    ζ = basis_change(layer.ζ(λ), inv(eul)) # HelperFunctions.jl

    M_mat(layer.d, ω, qₓ, ϵ, μ, ξ, ζ)
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

    if !allequal(diag(superstrate.ϵ(λ))) ||
       !allequal(diag(substrate.ϵ(λ)))   ||
       !allequal(diag(superstrate.μ(λ))) ||
       !allequal(diag(substrate.μ(λ)))
        throw(throw(ArgumentError("The superstrate and substrate need to be isotropic!")))
    end

    cos_αt, n_inc, n_out, qₓ = calculate_cos_αt(λ, α, strct) # HelperFunctions.jl

    K_sup = K_mat(n_inc, cos(α))
    K_sub = K_mat(n_out, cos_αt)

    if !isempty(layers)
	lay_Ms = calculate_layer_M.(λ, qₓ, layers)

	M = inv(K_sub) * prod(reverse(lay_Ms)) * K_sup
    else
	M = inv(K_sub) * K_sup
    end

    M
end
