"""
    nan_to_zero(x)

Return `0` if the input is `NaN`, othewise return `x`.
"""
nan_to_zero(x) = isnan(x) ? 0 : x

"""
    build_dir(f)

Creats path to files in build directory.
"""
build_dir(f) = joinpath(@__DIR__, "..", "deps", f)

##------------------------------------------------------------------------------
## Rotation matrix
##------------------------------------------------------------------------------


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


##------------------------------------------------------------------------------
## General basis change
##------------------------------------------------------------------------------

"""
   basis_change(M,T)

Change basis of matrix `M` using the basis transformation matrix `T`.
"""
basis_change(M, T) = inv(T) * M * T


##------------------------------------------------------------------------------
## Changing between linear and circular basis
##------------------------------------------------------------------------------


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
        throw(ArgumentError("$(repr(basis)) basis not defined. Only :linear and
:circular bases are defined."))
    end
end

##------------------------------------------------------------------------------
## Angle of transmission
##------------------------------------------------------------------------------


"""
    calculate_cos_αt(λ, α, strct)

Calculate the cosine outgoing angle for a layered structure `strct` and angle of
incidence `α`. Returns the tuple `(cos_αt, n_inc, n_out, qₓ)`, where `cos_αt` is
the cosine of the outgoing angle, `n_inc` is the refractive index of the
incident medium, `n_out` is the refractive index of the outgoing medium, and
`qₓ` is the normalized in-plane wavevector component.

We make sure to choose the correct branch when calculating the cosine for more
see Section 5 in [^1].

# Examples
```jldoctest
julia> S = LayeredStructure(superstrate=Layer(), substrate = SiC());
julia> calculate_cos_αt(12e-6, deg2rad(14), S)
(0.0014914047759539242 - 0.04676247666123382im, 1.0 + 0.0im, 0.1648892944698832
+ 5.166277293025589im, 0.24192189559966773 + 0.0im)
```

# References
[^1]: $(References["Byrnes"])
"""
function calculate_cos_αt(λ, α, strct)
    # Only works for isotropic in- and outgoing media
    @unpack superstrate, layers, substrate = strct

    n_inc = sqrt(Complex(superstrate.ϵ(λ)[1,1] * superstrate.μ(λ)[1,1]))
    n_out = sqrt(Complex(substrate.ϵ(λ)[1,1] * substrate.μ(λ)[1,1]))
    qₓ = n_inc * sin(α)

    # Naively calculate cos_αt
    sin_αt = (n_inc / n_out) * sin(α)
    cos_αt = sqrt(1 - sin_αt^2)

    # Make sure we choose the right branch
    ret = real(n_out * cos_αt)
    imt = imag(n_out * cos_αt)

    if real(n_out) > 0
        if imag(n_out) > 0
            if imt < 0
                cos_αt *= -1
            end
        elseif isapprox(imag(n_out),0,atol=1e-8)
            if isapprox(imt,0,atol=1e-14)
                if ret < 0
                    cos_αt *= -1
                end
            elseif imt < 0
                cos_αt *= -1
            end
        end
    # Backwards traveling wave
    elseif real(n_out) < 0 && imag(n_out) < 0
        if imt < 0
            cos_αt *= -1
	end
    end
    # This probably breaks for active media. The current approach would be to
    # pray we choose the right branch.

    cos_αt, n_inc, n_out, qₓ
end
