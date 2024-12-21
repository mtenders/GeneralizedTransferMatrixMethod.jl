"""
    struct Layer
        ϵ::Function
        d::Real
        θ::Real
        ϕ::Real
        ψ::Real
    end

Type that represents a layer in a layered structur. Boundary layers have d = 0.

### Fields

- `ϵ` -- Permitivity tensor in diagonal form (default: ϵ_vacuum).
- `d` -- Thickness of the layer `[m]` (default: 0).
- `θ` -- θ Euler angle `[rad]` (default: 0).
- `ϕ` -- ϕ Euler angle `[rad]` (default: 0).
- `ψ` -- ψ Euler angle `[rad]` (default: 0).

### Reference

The definition of the Euler matrix is taken from *Optical Waves in Layered Media
by Pochi Yeh*.
"""
@kwdef struct Layer
    ϵ::Function = ϵ_vacuum
    d::Real = 0
    θ::Real = 0
    ϕ::Real = 0
    ψ::Real = 0
end

"""
    struct LayeredStructure
        superstrate::Layer
        layers::AbstractVector{Layer}
        substrate::Layer
    end

Type that represents a layered structure.

### Fields

- superstrate -- Superstrate of the structure.
- layers      -- List of layers of the structure starting with the layer
                 underneath the superstrate and ending with the layer above the
                 substrate (default = `[]`).
- substrate   -- Subststrate of the structure.

"""
@kwdef struct LayeredStructure
    superstrate::Layer
    layers::AbstractVector{Layer} = []
    substrate::Layer
end

"""
    struct LayerProperties
        T::AbstractMatrix
        P::AbstractMatrix
        A::AbstractMatrix
        γ::AbstractMatrix
        q::AbstractVector
        Ψ::AbstractVector
        S::AbstractVector
        Δ::AbstractMatrix
        a::AbstractMatrix
        M::AbstractMatrix
    end

Type that holds all properties of a single layer in a structure. For boundary
layers P and T are identity matrices.

The properties of this type only make sense in combination with the wavelength λ
and the normalized in-plane wavevector component ζ.

### Fields

- T -- Transmission matrix of the layer.
- P -- Progation matrix of the layer.
- A -- Dynamical matrix of the layer.
- γ -- Normalized γ vectors in matrix form.
- q -- Four eigenvalues (z-components of the wavevector) of the layer.
- Ψ -- Four eigenmodes of the layer.
- S -- Poynting vectors of the four eigenmodes of the layer.
- Δ -- Δ matrix.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
@kwdef struct LayerProperties
    T::AbstractMatrix
    P::AbstractMatrix
    A::AbstractMatrix
    γ::AbstractMatrix
    q::AbstractVector
    Ψ::AbstractMatrix
    S::AbstractVector
    Δ::AbstractMatrix
    a::AbstractMatrix
    M::AbstractMatrix
end

"""
    struct StructureProperties
                Γ::AbstractMatrix
                superstrate::LayerProperties
                layers::AbstractVector{LayerProperties}
                substrate::LayerProperties
        end

Type that holds all properties of an entire layered structure.

The properties of this type only make sense in combination with the wavelength λ
and the normalized in-plane wavevector component ζ.

### Fields

- Γ           -- Full transfer matrix Γ* (see reference).
- superstrate -- Properties of the superstrate layer.
- layers      -- List of properties of the intermediate layers.
- substrate   -- Properties of the substrate layer.

### Reference

The definition is taken from [Passler and Paarmann
2017](https://doi.org/10.1364/JOSAB.34.002128) and [Passler and Paarmann 2019
(erratum)](https://doi.org/10.1364/JOSAB.36.003246).
"""
@kwdef struct StructureProperties
    Γ::AbstractMatrix
    superstrate::LayerProperties
    layers::AbstractVector{LayerProperties}
    substrate::LayerProperties
end

