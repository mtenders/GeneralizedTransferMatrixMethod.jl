##------------------------------------------------------------------------------
## HELPER FUNCTIONS
##------------------------------------------------------------------------------

"""
    nan_to_zero(x)

Return `0` if the input is `NaN`, othewise return `x`.
"""
nan_to_zero(x) = isnan(x) ? 0 : x

"""
   RvsI(x)

Check if the magnitude of the real part of a number `x` is bigger than the
maginute of the imaginary part. Return the bigger one of both. If they are the
same return the imaginary part.
"""
RvsI(x) = abs(real(x)) > abs(imag(x)) ? real(x) : imag(x)

##------------------------------------------------------------------------------
## EIGENVALUES
##------------------------------------------------------------------------------

"""
    eigen_trans_ref(Δ)

Return the four eigenvalues and vectors of the Δ matrix sorted by propagation
direction.

The first two entries represent the forward propagating (transmitted) waves, the
last two entries the backward propagating (reflected) waves.

### Input

- `Δ` -- Square matrix.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128).
"""
function eigen_trans_ref(Δ)
    by(λ) = (-imag(λ), -real(λ))
    q, Ψ = eigen(Δ, sortby=by)

    # The real part and the imaginary part should always have the same sign. To
    # avoid edge cases (materials with ϵ = diag(1000, -1000, 1000)), where the
    # real part is zero for one direction and the imaginary part is zero for the
    # other, we just use the bigger one for the sorting. The sorting within the
    # transmitted and reflected waves does not matter and will be taken care of
    # later.

    idx = sortperm(RvsI.(q), rev=true)
    q = q[idx]
    Ψ = Ψ[:,idx]

    return q, Ψ

    # # Check if the imaginary part is zero and sort again by real part to avoid
    # # sorting problems arising from +0im > -0im

    # # Within machine precision
    # if all(isapprox.(imag.(q),0.0,atol=eps(Float64), rtol=0))
    #     # Get sorted permutation, Sort by real part from biggest to smallest
    #     # number.
    #     idx = sortperm(real.(q), rev=true)
    #     q = q[idx]
    #     Ψ = Ψ[:,idx]
    # end
    # return q, Ψ
end


"""
    poynting_vector(Ψ, ζ, ϵ, μ)

Calculate the Poynting vector from the eigenvector Ψ of Δ.

### Input

- `Ψ` -- Eigenvector of Δ of the form (E_x, H_y, E_y, -H_x)ᵀ.
- `ζ` -- In-plane reduced wavevector kₓ/k₀ in the system.
- `ϵ` -- 3x3 permitivity tensor in the lab frame.
- `μ` -- Scalar permeability.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128).
"""
function poynting_vector(Ψ, ζ, ϵ, μ)
    M = M_mat(ϵ, μ)
    a = a_mat(M, ζ)
    E_x = Ψ[1]
    E_y = Ψ[3]
    H_x = -Ψ[4]
    H_y = Ψ[2]

    E_z = a[3,1] * E_x + a[3,2] * E_y + a[3,4] * H_x + a[3,5] * H_y
    H_z = a[6,1] * E_x + a[6,2] * E_y + a[6,4] * H_x + a[6,5] * H_y

    S = [E_y * H_z - E_z * H_y,
         E_z * H_x - E_x * H_z,
         E_x * H_y - E_y * H_x]
    return S
end


"""
    C_p(F₁, F₂)

Calculate the ordering parameter C from for two quantities F₁ and F₂.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128).
"""
C_p(F₁, F₂) = abs2(F₁) / (abs2(F₁) + abs2(F₂))


"""
    eigen_sorted(ζ, ϵ, μ, Δ)

Return the sorted eigenvalues, eigenvectors, and Poynting vectors of
Δ(ζ, ϵ, μ).

### Input

- `ζ` -- In-plane reduced wavevector kₓ/k₀ in the system.
- `ϵ` -- 3x3 permitivity tensor in the lab frame.
- `μ` -- Scalar permeability.
- `Δ` -- Δ matrix.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128).
"""
function eigen_sorted(ζ, ϵ, μ, Δ)
    # Sort into transmitted and reflected waves
    q, Ψ = eigen_trans_ref(Δ)
    # Calculate Poynting vectors for 4 Eigen vectors of Δ
    S = Vector{Vector{Complex}}(undef, 4) # Preallocate array for S
    C = Vector{Real}(undef,4) # Preallocate array for order parameter
    # using the Poynting vector
    for i = 1:4
        S[i] = poynting_vector(Ψ[:,i], ζ, ϵ, μ)
        C[i] = C_p(S[i][1], S[i][2])
    end

    # Non-birefringent: Override C using the definition involving electric
    # fields.
    # If the Poynting vector only has a component in the z direction C evaluates
    # to 0/0 = NaN.
    if C[1] ≈ C[2] || any(isnan, C)
        for i = 1:4
            C[i] = C_p(Ψ[1,i],Ψ[3,i])
        end
    end

    # Swap elements if necessary
    if C[2] > C[1]
        q[1:2] = reverse(q[1:2])
        Ψ[:,1:2] = reverse(Ψ[:,1:2], dims=2)
        S[1:2] = reverse(S[1:2])
    end
    if C[4] > C[3]
        q[3:4] = reverse(q[3:4])
        Ψ[:,3:4] = reverse(Ψ[:,3:4], dims=2)
        S[3:4] = reverse(S[3:4])
    end
    return q, Ψ, S
end


##------------------------------------------------------------------------------
## PROPERTIES
##------------------------------------------------------------------------------


"""
    calculate_layer_properties(layer::Layer, ζ::Real, λ::Real)

Calculate the properties of a single layer in a layered structure and returns a
struct of type LayerProperties.

### Input

- `layer` -- Layer.
- `ζ`     -- In-plane reduced wavevector kₓ/k₀ in the system.
- `λ`     -- Wavelength `[m]`.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
function calculate_layer_properties(layer::Layer, ζ::Real, λ::Real)
    μ = 1.0
    # ϵ's in lab frame intermediate layers
    eul = euler_mat(layer)
    ϵ_lay = inv(eul) * layer.ϵ(λ) * eul # TODO check again
    # M matrix
    M = M_mat(ϵ_lay, μ)
    # a matrix
    a = a_mat(M, ζ)
    # Δ matrix
    Δ = Δ_mat(ζ, ϵ_lay, μ, M, a)
    # Four sorted eigen values, eigen vectors, and Poynting vectors.
    q, Ψ, S = eigen_sorted(ζ, ϵ_lay, μ, Δ)
    # γ matrix
    γ = γ_mat(q, ζ, ϵ_lay, μ)
    # Dynamical matrix
    A = dynamical_mat(q, ζ, ϵ_lay, μ, γ)
    # Propagation matrix
    P = propagation_mat(q, λ, layer.d)
    # Transmission matrix
    T = A * P * inv(A)
    # Fill structure
    return LayerProperties(T, P, A, γ, q, Ψ, S, Δ, a, M)
end


"""
    calculate_structure_properties(ζ::Real, λ::Real, strct::LayeredStructure)

Calculate the properties of a layered structure and returns a struct of type
`StructureProperties`.

### Input

- `ζ`     -- In-plane reduced wavevector kₓ/k₀ in the system.
- `λ`     -- Wavelength `[m]`.
- `strct` -- Layered structure.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
function calculate_structure_properties(ζ::Real, λ::Real, strct::LayeredStructure)
    @unpack superstrate, layers, substrate = strct

    sup_props = calculate_layer_properties(superstrate, ζ, λ)
    sub_props = calculate_layer_properties(substrate, ζ, λ)

    if !isempty(layers)
        lay_props = calculate_layer_properties.(layers, ζ, λ)
        # Asup⁻¹ * ∏Tᵢ * Asub
        Γ = inv(sup_props.A) * prod(getproperty.(lay_props, :T)) * sub_props.A
    else
        lay_props = []
        Γ = inv(sup_props.A) * sub_props.A
    end

    Γstar = inv(Λ) * Γ * Λ

    return StructureProperties(Γstar, sup_props, lay_props, sub_props)
end


"""
    calculate_structure_properties(ζ::Real, λ::Unitful.Length,
                                   strct::LayeredStructure)

Calculate the properties of a layered structure using Unitful quantities and
returns a struct of type `StructureProperties`.

### Input

- `ζ`     -- In-plane reduced wavevector kₓ/k₀ in the system.
- `λ`     -- Wavelength.
- `strct` -- Layered structure.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
function calculate_structure_properties(ζ::Real, λ::Unitful.Length,
                                        strct::LayeredStructure)
    calculate_structure_properties(ζ, ustrip(u"m",λ), strct)
end


##------------------------------------------------------------------------------
## GLOBAL PROPERTIES
##------------------------------------------------------------------------------


"""
    reflection_coeffs(Γ)

Calculate the reflection coefficients from the full transfer matrix Γ (Γ* in the
reference). Returns the reflection coefficients in the order rₚₚ, rₛₛ, rₚₛ,
rₛₚ.

### Input

- `Γ` -- Full transfer matrix in Yeh's formalism.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
function reflection_coeffs(Γ)
    denom1 = Γ[1,1] - Γ[1,3] * Γ[3,1] / Γ[3,3]
    denom2 = Γ[3,3] - Γ[1,3] * Γ[3,1] / Γ[1,1]

    r_pp = (Γ[2,1] - Γ[2,3] * Γ[3,1] / Γ[3,3]) / denom1
    r_ss = (Γ[4,3] - Γ[4,1] * Γ[1,3] / Γ[1,1]) / denom2
    r_ps = (Γ[4,1] - Γ[4,3] * Γ[3,1] / Γ[3,3]) / denom1
    r_sp = (Γ[2,3] - Γ[2,1] * Γ[1,3] / Γ[1,1]) / denom2

    return r_pp, r_ss, r_ps, r_sp
end
"""
    reflection_coeffs(props::StructureProperties)

Calculate the reflection coefficients from a `StructureProperties` type. Returns the
reflection coefficients for different polarisations in the order rₚₚ, rₛₛ, rₚₛ, rₛₚ.

### Input

- `props` -- Layered material reperesented by a `StructureProperties` type.

"""
function reflection_coeffs(props::StructureProperties)
    reflection_coeffs(props.Γ)
end

"""
    transmission_coeffs(Γ)

Calculate the transmission coefficients from the full transfer matrix Γ (Γ* in
the reference). Returns the reflection coefficients in the order tₚₚ, tₛₛ, tₚₛ,
tₛₚ.

### Input

- `Γ` -- Full transfer matrix in Yeh's formalism.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
function transmission_coeffs(Γ)
    denom = Γ[3,3] * Γ[1,1] - Γ[1,3] * Γ[3,1]

    t_pp = 1 / (Γ[1,1] - Γ[1,3] * Γ[3,1] / Γ[3,3])
    t_ss = 1 / (Γ[3,3] - Γ[1,3] * Γ[3,1] / Γ[1,1])

    # TODO
    # Filter case where denom is NaN but Γ[3,1] is 0
    t_ps = iszero(Γ[3,1]) ? 0 : -Γ[3,1] / denom
    t_sp = iszero(Γ[1,3]) ? 0 : -Γ[1,3] / denom

    return t_pp, t_ss, t_ps, t_sp
end
"""
    transmission_coeffs(props::StructureProperties)

Calculate the transmission coefficients from a `StructureProperties` type. Returns the
transmission coefficients for different polarisations in the order tₚₚ, tₛₛ, tₚₛ, tₛₚ.

### Input

- `props` -- Layered material reperesented by a `StructureProperties` type.

"""
function transmission_coeffs(props::StructureProperties)
    transmission_coeffs(props.Γ)
end

"""
    reflection(r)

Calculate the reflection from the reflection coefficient. (Just square moduli)

### Input

- `r` -- Reflection coefficient.

"""
function reflection(r)
    abs2(r)
end
"""
    reflection(props::StructureProperties)

Calculate the reflection from a `StructureProperties` type. Returns the
reflection for different polarisations in the order Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ.

### Input

- `props` -- Layered material reperesented by a `StructureProperties` type.

"""
function reflection(props::StructureProperties)
    reflection.(reflection_coeffs(props))
end


"""
    transmission(ζ, props::StructureProperties)

Calculate the total transmission from a `StructureProperties` type. Returns
the transmission in the order Tₚ, Tₛ, where the indices indicate the incoming
polarisation.

### Input

- `ζ`     -- In-plane reduced wavevector kₓ/k₀ of the system.
- `props` -- Layered material reperesented by a `StructureProperties` type.

"""
function transmission(ζ, props::StructureProperties)
    # INCOMING FIELDS
    # Save 4 k-vectors in a matrix (superstrate)
    k_sup = [ ones(4)*ζ zeros(4) props.superstrate.q ]

    E_inc_pin = props.superstrate.γ[1,:]
    E_inc_sin = props.superstrate.γ[2,:]

    S_inc_pin = 0.5 * real(E_inc_pin × conj(k_sup[1,:] × E_inc_pin))
    S_inc_sin = 0.5 * real(E_inc_sin × conj(k_sup[2,:] × E_inc_sin))

    # OUTGOING FIELDS
    # Save 4 k-vectors in a matrix (substrate)
    k_sub = [ ones(4)*ζ zeros(4) props.substrate.q]
    t_pp, t_ss, t_ps, t_sp = transmission_coeffs(props.Γ)

    E_out_pin = t_pp * props.substrate.γ[1,:] + t_ps * props.substrate.γ[2,:] #p-in, p or s out
    E_out_sin = t_sp * props.substrate.γ[1,:] + t_ss * props.substrate.γ[2,:] #s-in, p or s out

    S_out_pin = 0.5 * real(E_out_pin × conj(k_sub[1,:] × E_out_pin))
    S_out_sin = 0.5 * real(E_out_sin × conj(k_sub[2,:] × E_out_sin))

    # Transmission from the z-component of the Poynting vectors
    T_p = S_out_pin[3]/S_inc_pin[3]
    T_s = S_out_sin[3]/S_inc_sin[3]

    return T_p, T_s
end
