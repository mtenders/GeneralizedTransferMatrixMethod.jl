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
    [ ϵ          zeros(3,3);
      zeros(3,3) μ * I]
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
    a = zeros(Complex, 6,6)
    b = M[3,3] * M[6,6] - M[3,6] * M[6,3]

    a[3, 1] = (M[6,1] * M[3,6] - M[3,1] * M[6,6]) / b
    a[3, 2] = ((M[6,2] - ζ) * M[3,6] - M[3,2] * M[6,6]) / b
    a[3, 4] = (M[6,4] * M[3,6] - M[3,4] * M[6,6]) / b
    a[3, 5] = (M[6,5] * M[3,6] - (M[3,5] + ζ) * M[6,6]) / b
    a[6, 1] = (M[6,3] * M[3,1] - M[3,3] * M[6,1]) / b
    a[6, 2] = (M[6,3] * M[3,2] - M[3,3] * (M[6,2] - ζ)) / b
    a[6, 4] = (M[6,3] * M[3,4] - M[3,3] * M[6,4]) / b
    a[6, 5] = (M[6,3] * (M[3,4] + ζ) - M[3,3] * M[6,5]) / b
    return a
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
    Δ = Matrix{Complex}(undef, 4,4)

    Δ[1,1] = M[5,1] + (M[5,3] + ζ) * a[3,1] + M[5,6] * a[6,1]
    Δ[1,2] = M[5,5] + (M[5,3] + ζ) * a[3,5] + M[5,6] * a[6,5]
    Δ[1,3] = M[5,2] + (M[5,3] + ζ) * a[3,2] + M[5,6] * a[6,2]
    Δ[1,4] = -M[5,4] - (M[5,3] + ζ) * a[3,4] - M[5,6] * a[6,4]
    Δ[2,1] = M[1,1] + M[1,3] * a[3,1] + M[1,6] * a[6,1]
    Δ[2,2] = M[1,5] + M[1,3] * a[3,5] + M[1,6] * a[6,5]
    Δ[2,3] = M[1,2] + M[1,3] * a[3,2] + M[1,6] * a[6,2]
    Δ[2,4] = -M[1,4] - M[1,3] * a[3,4] + M[1,6] * a[6,4]
    Δ[3,1] = -M[1,4] - M[4,3] * a[3,1] - M[4,6] * a[6,1]
    Δ[3,2] = -M[4,5] - M[4,3] * a[4,3] - M[4,6] * a[6,5]
    Δ[3,3] = -M[4,2] - M[4,3] * a[3,2] - M[4,6] * a[6,2]
    Δ[3,4] = M[4,4] + M[4,3] * a[3,4] + M[4,6] * a[6,4]
    Δ[4,1] = M[2,1] + M[2,3] * a[3,1] + (M[2,6] - ζ) * a[6,1]
    Δ[4,2] = M[2,5] + M[2,3] * a[3,5] + (M[2,6] - ζ) * a[6,5]
    Δ[4,3] = M[2,2] + M[2,3] * a[3,2] + (M[2,6] - ζ) * a[6,2]
    Δ[4,4] = -M[2,4] - M[2,3] * a[3,4] + (M[2,6] - ζ) * a[6,4]
    return Δ
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
    γ = Matrix{Complex}(undef, 4, 3)
    γ[1,1] = γ[2,2] = γ[4,2] = 1 + 0im
    γ[3,1] = -1 + 0im

    if q[1] ≈ q[2]
        γ[1,2] = 0 + 0im
        γ[1,3] = -(μ * ϵ[3,1] + ζ * q[1]) / (μ * ϵ[3,3] - ζ^2)
        γ[2,1] = 0 + 0im
        γ[2,3] = -μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
    else # q[1] ≉ q[2]
        γ12_num = μ * ϵ[2,3] * (μ * ϵ[3,1] + ζ * q[1]) - μ * ϵ[2,1] * (μ * ϵ[3,3] - ζ^2)
        γ12_denom = (μ * ϵ[3,3] - ζ^2) * (μ * ϵ[2,2] - ζ^2 - q[1]^2) - μ^2 * ϵ[2,3] * ϵ[3,2]
        γ[1,2] = γ12_num / γ12_denom
        γ13_a = -(μ * ϵ[3,1] + ζ * q[1]) / (μ * ϵ[3,3] - ζ^2)
        γ13_b = -μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
        γ[1,3] = γ13_a + γ13_b * γ[1,2]
        γ21_num = μ * ϵ[3,2] * (μ * ϵ[1,3] + ζ * q[2]) - μ * ϵ[1,2] * (μ * ϵ[3,3] - ζ^2)
        γ21_denom = (μ * ϵ[3,3] - ζ^2) * (μ * ϵ[1,1] - q[2]^2) - (μ * ϵ[1,3] + ζ * q[2]) * (μ * ϵ[3,1] + ζ * q[2])
        γ[2,1] = γ21_num / γ21_denom
        γ23_a = -(μ * ϵ[3,1] + ζ * q[2]) / (μ * ϵ[3,3] - ζ^2)
        γ23_b = -μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
        γ[2,3] = γ23_a * γ[2,1] + γ23_b
    end

    if q[3] ≈ q[4]
        γ[3,2] = 0 + 0im
        γ[3,3] = (μ * ϵ[3,1] + ζ * q[3]) / (μ * ϵ[3,3] - ζ^2)
        γ[4,1] = 0 + 0im
        γ[4,3] = -μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
    else # q[3] ≉ q[4]
        # in the paper it's written wrong as μ * ϵ[3,3] + ζ^2, the Xu paper uses μ * ϵ[3,3] - ζ^2
        γ32_num = μ * ϵ[2,1] * (μ * ϵ[3,3] - ζ^2) - μ * ϵ[2,3] * (μ * ϵ[3,1] + ζ * q[3])
        γ32_denom = (μ * ϵ[3,3] - ζ^2) * (μ * ϵ[2,2] - ζ^2 - q[3]^2) - μ^2 * ϵ[2,3] * ϵ[3,2]
        γ[3,2] = γ32_num / γ32_denom
        γ33_a = (μ * ϵ[3,1] + ζ * q[3]) / (μ * ϵ[3,3] - ζ^2)
        γ33_b = μ * ϵ[3,2] / (μ * ϵ[3,3] - ζ^2)
        γ[3,3] = γ33_a + γ33_b * γ[3,2]
        γ41_num = μ * ϵ[3,2] * (μ * ϵ[1,3] + ζ * q[4]) - μ * ϵ[1,2] * (μ * ϵ[3,3] - ζ^2)
        γ41_denom = (μ * ϵ[3,3] - ζ^2) * (μ * ϵ[1,1] - q[4]^2) - (μ * ϵ[1,3] + ζ * q[4]) * (μ * ϵ[3,1] + ζ * q[4])
        γ[4,1] = γ41_num / γ41_denom
        γ43_num = -(μ * ϵ[3,1] + ζ * q[4]) * γ[4,1] - μ * ϵ[3,2]
        γ43_denom =  μ * ϵ[3,3] - ζ^2
        γ[4,3] = γ43_num / γ43_denom
    end

    # Normalize the γ vectors
    for i = 1:4
        γ[i,:] = normalize(γ[i,:])
    end
    return γ
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
    A = Matrix{Complex}(undef, 4, 4)
    A[1,:] = γ[:,1]
    A[2,:] = γ[:,2]
    A[3,:] = (q .* γ[:,1] .- ζ .* γ[:,3]) / μ
    A[4,:] = q .* γ[:,2] / μ
    return SMatrix{4,4}(A)
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
    SMatrix{4,4}(Diagonal(exp.(-1im * 2π * q * d / λ)))
end
