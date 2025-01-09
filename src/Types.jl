"""
    struct Layer
        ϵ::Function = ϵ_vacuum
        μ::Function = μ_vacuum
        ξ::Function = ξ_vacuum
        ζ::Function = ζ_vacuum
        d::Real = 0
        θ::Real = 0
        ϕ::Real = 0
        ψ::Real = 0
    end

Type that represents a layer in a layered structur. Boundary layers have `d = 0`.

# Arguments

- `ϵ`: Relative permitivity tensor (default: `ϵ_vacuum`).
- `μ`: Relative permeability tensor (default: `μ_vacuum`).
- `ξ`, `ζ`: Optical rotation tensors (default: `ξ_vacuum`, `ζ_vacuum`).
- `d`: Thickness of the layer `[m]` (default: 0).
- `θ`: θ Euler angle `[rad]` (default: 0).
- `ϕ`: ϕ Euler angle `[rad]` (default: 0).
- `ψ`: ψ Euler angle `[rad]` (default: 0).

The definition of the Euler angles is taken from [^1].

# References
[^1]: $(References["Yeh"])
"""
@kwdef struct Layer
    ϵ::Function = ϵ_vacuum
    μ::Function = μ_vacuum
    ξ::Function = ξ_vacuum
    ζ::Function = ζ_vacuum
    d::Real = 0
    θ::Real = 0
    ϕ::Real = 0
    ψ::Real = 0
end
Base.broadcastable(f::Layer) = Ref(f)

"""
    struct LayeredStructure
        superstrate::Layer
        layers::AbstractVector{Layer} = []
        substrate::Layer
    end

Type that represents a layered structure.

# Arguments

- superstrate: Superstrate of the structure.
- layers: List of layers of the structure starting with the layer underneath the
          superstrate and ending with the layer above the substrate (default = `[]`).
- substrate: Subststrate of the structure.

"""
@kwdef struct LayeredStructure
    superstrate::Layer
    layers::AbstractVector{Layer} = []
    substrate::Layer
end
Base.broadcastable(f::LayeredStructure) = Ref(f)
