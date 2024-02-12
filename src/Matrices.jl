"""
    Λ

Λ matrix (see reference). Changes the order of the Γ matrix to follow Yeh's
formalism.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128).
"""
const Λ = @SMatrix [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]


"""
    euler_mat(θ, ϕ, ψ)

Calculate the Euler matrix from the Euler angles `θ`, `ϕ`, and `ψ`.
"""
function euler_mat(θ, ϕ, ψ)
    @SMatrix [cos(ψ)*cos(ϕ)-cos(θ)*sin(ϕ)*sin(ψ) -sin(ψ)*cos(ϕ)-cos(θ)*sin(ϕ)*cos(ψ)  sin(θ)*sin(ϕ);
              cos(ψ)*sin(ϕ)+cos(θ)*cos(ϕ)*sin(ψ) -sin(ψ)*sin(ϕ)+cos(θ)*cos(ϕ)*cos(ψ) -sin(θ)*cos(ϕ);
              sin(θ)*sin(ψ)                       sin(θ)*cos(ψ)                       cos(θ)]
end
"""
        euler_mat(layer::Layer)

Calculate the Euler matrix for a Layer.
"""
function euler_mat(layer::Layer)
    @unpack θ, ϕ, ψ = layer
    euler_mat(θ, ϕ, ψ)
end


"""
    M_mat(ϵ, μ)

Calculate the M matrix (see reference).

### Input

- `ϵ` -- 3x3 permitivity tensor in the lab frame.
- `μ` -- Scalar permeability.

### Notes

Only non-optically active media with isotropic permeability are considered:
```math
\\bar\\mu = \\mu\\bar1 \\quad
\\bar\\rho_1 = \\bar\\rho_2 = \\bar0
```

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128).
"""
function M_mat(ϵ, μ)
    # I didn't know how to construct this matrix in a
    # more elegant way using @SMatrix
    @SMatrix [
	ϵ[1,1] ϵ[1,2] ϵ[1,3] 0 0 0
	ϵ[2,1] ϵ[2,2] ϵ[2,3] 0 0 0
	ϵ[3,1] ϵ[3,2] ϵ[3,3] 0 0 0
	0 0 0 μ 0 0
	0 0 0 0 μ 0
	0 0 0 0 0 μ
    ]
end


"""
    a_mat(M, ζ)

Calculate the a matrix (see reference).

### Input

- `M` -- 6x6 M matrix (see reference).
- `ζ` -- In-plane reduced wavevector kₓ/k₀ in the system.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128).
"""
function a_mat(M, ζ)
    b = M[3,3] * M[6,6] - M[3,6] * M[6,3]

    a31 = (M[6,1] * M[3,6] - M[3,1] * M[6,6]) / b
    a32 = ((M[6,2] - ζ) * M[3,6] - M[3,2] * M[6,6]) / b
    a34 = (M[6,4] * M[3,6] - M[3,4] * M[6,6]) / b
    a35 = (M[6,5] * M[3,6] - (M[3,5] + ζ) * M[6,6]) / b
    a61 = (M[6,3] * M[3,1] - M[3,3] * M[6,1]) / b
    a62 = (M[6,3] * M[3,2] - M[3,3] * (M[6,2] - ζ)) / b
    a64 = (M[6,3] * M[3,4] - M[3,3] * M[6,4]) / b
    a65 = (M[6,3] * (M[3,4] + ζ) - M[3,3] * M[6,5]) / b

    @SMatrix [
	0 0 0 0 0 0
	0 0 0 0 0 0
	a31 a32 0 a34 a35 0
	0 0 0 0 0 0
	0 0 0 0 0 0
	a61 a62 0 a64 a65 0
    ]
end


"""
    Δ_mat(ζ, ϵ, μ, M, a)

Calculate the Δ matrix (see reference).

### Input

- `ζ` -- In-plane reduced wavevector kₓ/k₀ in the system.
- `ϵ` -- 3x3 permitivity tensor in the lab frame.
- `μ` -- Scalar permeability.
- `M` -- M matrix.
- `a` -- a matrix.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128).
"""
function Δ_mat(ζ, ϵ, μ, M, a)
    Δ11 = M[5,1] + (M[5,3] + ζ) * a[3,1] + M[5,6] * a[6,1]
    Δ12 = M[5,5] + (M[5,3] + ζ) * a[3,5] + M[5,6] * a[6,5]
    Δ13 = M[5,2] + (M[5,3] + ζ) * a[3,2] + M[5,6] * a[6,2]
    Δ14 = -M[5,4] - (M[5,3] + ζ) * a[3,4] - M[5,6] * a[6,4]
    Δ21 = M[1,1] + M[1,3] * a[3,1] + M[1,6] * a[6,1]
    Δ22 = M[1,5] + M[1,3] * a[3,5] + M[1,6] * a[6,5]
    Δ23 = M[1,2] + M[1,3] * a[3,2] + M[1,6] * a[6,2]
    Δ24 = -M[1,4] - M[1,3] * a[3,4] + M[1,6] * a[6,4]
    Δ31 = -M[1,4] - M[4,3] * a[3,1] - M[4,6] * a[6,1]
    Δ32 = -M[4,5] - M[4,3] * a[4,3] - M[4,6] * a[6,5]
    Δ33 = -M[4,2] - M[4,3] * a[3,2] - M[4,6] * a[6,2]
    Δ34 = M[4,4] + M[4,3] * a[3,4] + M[4,6] * a[6,4]
    Δ41 = M[2,1] + M[2,3] * a[3,1] + (M[2,6] - ζ) * a[6,1]
    Δ42 = M[2,5] + M[2,3] * a[3,5] + (M[2,6] - ζ) * a[6,5]
    Δ43 = M[2,2] + M[2,3] * a[3,2] + (M[2,6] - ζ) * a[6,2]
    Δ44 = -M[2,4] - M[2,3] * a[3,4] + (M[2,6] - ζ) * a[6,4]

    @SMatrix[
	Δ11 Δ12 Δ13 Δ14
	Δ21 Δ22 Δ23 Δ24
	Δ31 Δ32 Δ33 Δ34
	Δ41 Δ42 Δ43 Δ44
    ]
end


"""
    γ_mat(q, ζ, ϵ, μ)

Calculate the normalized γ vectors and return them in matrix form.

### Input

- `q` -- Sorted eigenvalues of Δ.
- `ζ` -- In-plane reduced wavevector kₓ/k₀ in the system.
- `ϵ` -- 3x3 permitivity tensor in the lab frame.
- `μ` -- Scalar permeability.

### Notes

The field vectors are given by the rows of the matrix, so the indexing scheme
follows the one in the references.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
function γ_mat(q, ζ, ϵ, μ)
    γ11 = γ22 = γ42 = 1 + 0im
    γ31 = -1 + 0im

    if q[1] ≈ q[2]
        γ12 = 0 + 0im
        γ13 = -(μ * ϵ[3,1] + ζ * q[1]) / (μ * ϵ[3,3] - ζ^2)
        γ21 = 0 + 0im
        γ23 = -μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
    else # q[1] ≉ q[2]
        γ12_num = μ * ϵ[2,3] * (μ * ϵ[3,1] + ζ * q[1]) - μ * ϵ[2,1] * (μ * ϵ[3,3] - ζ^2)
        γ12_denom = (μ * ϵ[3,3] - ζ^2) * (μ * ϵ[2,2] - ζ^2 - q[1]^2) - μ^2 * ϵ[2,3] * ϵ[3,2]
        γ12 = γ12_num / γ12_denom
        γ13_a = -(μ * ϵ[3,1] + ζ * q[1]) / (μ * ϵ[3,3] - ζ^2)
        γ13_b = -μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
        γ13 = γ13_a + γ13_b * γ12
        γ21_num = μ * ϵ[3,2] * (μ * ϵ[1,3] + ζ * q[2]) - μ * ϵ[1,2] * (μ * ϵ[3,3] - ζ^2)
        γ21_denom = (μ * ϵ[3,3] - ζ^2) * (μ * ϵ[1,1] - q[2]^2) - (μ * ϵ[1,3] + ζ * q[2]) * (μ * ϵ[3,1] + ζ * q[2])
        γ21 = γ21_num / γ21_denom
        γ23_a = -(μ * ϵ[3,1] + ζ * q[2]) / (μ * ϵ[3,3] - ζ^2)
        γ23_b = -μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
        γ23 = γ23_a * γ21 + γ23_b
    end

    if q[3] ≈ q[4]
        γ32 = 0 + 0im
        γ33 = (μ * ϵ[3,1] + ζ * q[3]) / (μ * ϵ[3,3] - ζ^2)
        γ41 = 0 + 0im
        γ43 = -μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
    else # q[3] ≉ q[4]
        # in the paper it's written wrong as μ * ϵ[3,3] + ζ^2, the Xu paper uses μ * ϵ[3,3] - ζ^2
        γ32_num = μ * ϵ[2,1] * (μ * ϵ[3,3] - ζ^2) - μ * ϵ[2,3] * (μ * ϵ[3,1] + ζ * q[3])
        γ32_denom = (μ * ϵ[3,3] - ζ^2) * (μ * ϵ[2,2] - ζ^2 - q[3]^2) - μ^2 * ϵ[2,3] * ϵ[3,2]
        γ32 = γ32_num / γ32_denom
        γ33_a = (μ * ϵ[3,1] + ζ * q[3]) / (μ * ϵ[3,3] - ζ^2)
        γ33_b = μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
        γ33 = γ33_a + γ33_b * γ32
        γ41_num = μ * ϵ[3,2] * (μ * ϵ[1,3] + ζ * q[4]) - μ * ϵ[1,2] * (μ * ϵ[3,3] - ζ^2)
        γ41_denom = (μ * ϵ[3,3] - ζ^2) * (μ * ϵ[1,1] - q[4]^2) - (μ * ϵ[1,3] + ζ * q[4]) * (μ * ϵ[3,1] + ζ * q[4])
        γ41 = γ41_num / γ41_denom
        γ43_num = -(μ * ϵ[3,1] + ζ * q[4]) * γ41 - μ * ϵ[3,2]
        γ43_denom =  μ * ϵ[3,3] - ζ^2
        γ43 = γ43_num / γ43_denom
    end

    # I do the normalization by hand so I can use @SMatrix. There is probably a more elegant way to do this.
    γ1norm = norm([γ11, γ12, γ13])
    γ2norm = norm([γ21, γ22, γ23])
    γ3norm = norm([γ31, γ32, γ33])
    γ4norm = norm([γ41, γ42, γ43])
    @SMatrix [
	γ11/γ1norm γ12/γ1norm γ13/γ1norm
	γ21/γ2norm γ22/γ2norm γ23/γ2norm
	γ31/γ3norm γ32/γ3norm γ33/γ3norm
	γ41/γ4norm γ42/γ4norm γ43/γ4norm
    ]
end


"""
    dynamical_mat(q, ζ, ϵ, μ γ)

Calculate the dynamical matrix A of an interface. Returns an array of type
`SMatrix` for faster matrix operations.

### Input

- `q` -- Sorted eigenvalues of Δ.
- `ζ` -- In-plane reduced wavevector kₓ/k₀ in the system.
- `ϵ` -- 3x3 permitivity tensor in the lab frame.
- `μ` -- Scalar permeability.
- `γ` -- γ Matrix.


### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
function dynamical_mat(q, ζ, ϵ, μ, γ)
    A1 = transpose(γ[:,1])
    A2 = transpose(γ[:,2])
    A3 = transpose((q .* γ[:,1] .- ζ .* γ[:,3]) / μ)
    A4 = transpose(q .* γ[:,2] / μ)

    [A1; A2;  A3; A4]
end


"""
    propagation_mat(q, λ, d)

Calculate the propagation matrix P. Returns an array of type `SMatrix` for
faster matrix operations.

### Input

- `q` -- Sorted eigenvalues of Δ.
- `λ` -- Wavelength `[m]`.
- `d` -- Thickness of the material `[m]`.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
function propagation_mat(q, λ, d)
    @SMatrix [
	exp.(-1im * 2π * q[1] * d / λ) 0 0 0
	0 exp.(-1im * 2π * q[2] * d / λ) 0 0
	0 0 exp.(-1im * 2π * q[3] * d / λ) 0
	0 0 0 exp.(-1im * 2π * q[4] * d / λ)
    ]
end
