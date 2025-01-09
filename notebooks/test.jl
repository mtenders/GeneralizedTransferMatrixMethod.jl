### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 686524a8-d6ed-4ee9-9f27-fb1502b613d1
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate()
    using PlutoLinks: @revise 
	using LinearAlgebra
end

# ╔═╡ 68f96782-42d1-470f-bf10-6eefcd77186c
using Unitful: μm, m, °

# ╔═╡ 70e0c156-ce29-4927-94c7-64769c0261a0
using Unzip

# ╔═╡ 5e39aff5-7e98-4f68-aed5-a81ace47714f
using Plots

# ╔═╡ 1f69912d-4737-46f9-878a-562d02828ebb
@revise using GeneralizedTransferMatrixMethod

# ╔═╡ 107c43c8-6955-46b9-9849-4d45ad308535
md"""
# v0.2: The MacKay submodule

**ToDos**
- ✓ replace Parameters.jl's `@with_kw` with julia internal `@kw_def`
- ✓ implement MacKay code
  - ✓ Examples in Docstrings (``jldoctest`` macro?)
    - ✓ https://docs.julialang.org/en/v1/manual/documentation/
- ✓ Benchmark speed (new code is faster)
- Benchmark MacKay against old code
  - ✓ Waveplate paper
  - ✓ Chirality paper
  - ✓ Michela's paper (already exists)
  - ✓ Passler Paarman Otto geometry
- ✓ circular basis
- change basis change function so you dont have to insert inv(T)
- check the euler rotation against vera's code
- check what to export carefully (e.g. reflection_coeffs(M) necessary?)
- Rework Documentation (mention that it's only for isotropic in and out-going media)
  - Mention that circular basis is different from MacKay
  - make sure to correctly define r_ij and so on
- reexport still needed?
- prepare message that informs people about breaking changes upon package update (see Pluto code for how to do it)
- Unitful Extension
  - How does this work with the @permittivity macro
  - Test if Unitful is loaded
  - Define Layer for Unitful maybe then it's not needed inside the macro, documentation with @doc macro?
  - `Layer()` should also work with units
  - use k in 1/cm ?
- Makie Extension
- Plots Extension
- Rewrite documentation
- Include examples from papers
  - Write into Pluto notebooks (or maybe generate exmaples in documentation from notebooks)
    - Custom script that starts pluto with custom startpage?
- Include default permittivity models
- Put Readme in docstring
- Readme for src folder
- Typing??
- permeability and optical rotation tensors are included but NOT tested
- redo serious benchmarking

**Future features**
- `calculate_q_z`, eigenvalue of M
- dispersion plot
- Implement jump condition mention in MacKay p. 38
- Incoherent layers (Byrnes)
- Absorption of single layers
- Calculate electric field in structure
- `reverse` for LayeredStructure

**Tests**
- Reflection at 45° incidence is very commonly used for making 90° turns. For the case of light traversing from a less dense medium into a denser one at 45° incidence (θ = 45°), it follows algebraically from the above equations that Rp equals the square of Rs: R p = R s 2 ``{\displaystyle R_{\text{p}}=R_{\text{s}}^{2}}``
- Make sure the calculation of the outgoing angle is correct. Worst case we restrict it to lossless incoming and outgoing media (see Byrnes)


**Note for later**
- We transform t and r the same way. MacKay adds a minus
"""

# ╔═╡ 92640359-3f64-4b9c-ac7c-97aa58cd6478
TMM = GeneralizedTransferMatrixMethod

# ╔═╡ 984350af-8423-4238-8742-9e5a78293693
n_Si = 3.42

# ╔═╡ 214f6ace-39d3-4ad7-a05b-415fef8ca763
@permittivity "Si" λ -> n_Si^2 * Diagonal(ones(3));

# ╔═╡ 69731a15-62c7-4a93-b183-ce94f2132a37
Air = Layer()

# ╔═╡ 5f2a7902-8cfa-45a0-8e36-172fd3542ba6
md"""
## Waveplate paper
"""

# ╔═╡ 94acbf39-b4d7-47df-b121-95292a75b3a5
md"""
### Setup
"""

# ╔═╡ a0d2f37c-b411-4ee9-9913-e96204e0df64
function ϵ_x_MoO₃(λ)
    # Convert λ in meter to frequency in cm⁻¹
    f = 1 / (λ * 1e2)
    # Parameters from paper
    ϵ∞_x = 4.40
    fₗₒ_x₁, fₗₒ_x₂ = 606.9, 996.7 # [cm⁻¹]
    fₜₒ_x₁, fₜₒ_x₂ = 589.5, 811.5 # [cm⁻¹]
    γ_x₁, γ_x₂ = 51.53, 17.5# [cm⁻¹]

    return (ϵ∞_x * TMM.lorentz_osc(f, fₗₒ_x₁, fₜₒ_x₁, γ_x₁)
            * TMM.lorentz_osc(f, fₗₒ_x₂, fₜₒ_x₂, γ_x₂))
end

# ╔═╡ 78a20169-1ae5-4175-b139-9e11adf8f6d9
function ϵ_y_MoO₃(λ)
    # Convert λ in meter to frequency in cm⁻¹
    f = 1 / (λ * 1e2)
    # Parameters from paper
    ϵ∞_y = 4.86#4.87
    fₗₒ_y = 888.5#888.5 # [cm⁻¹]
    fₜₒ_y = 516.3#516.3 # [cm⁻¹]
    γ_y = 18.06#18.07 # [cm⁻¹]

    return (ϵ∞_y * TMM.lorentz_osc(f, fₗₒ_y, fₜₒ_y, γ_y))
end

# ╔═╡ 8193e9a3-1544-4bbd-a52a-7d4d7e47b769
@permittivity "MoO₃" λ -> Diagonal([ϵ_x_MoO₃(λ), ϵ_y_MoO₃(λ), ϵ_x_MoO₃(λ)]);

# ╔═╡ 94464fa0-21b1-4e4b-bfcb-98158d53e0a0
"""
	TStructure(thickness, ϕ=45°)

Semi-infinite Silicon as incident medium with α-MoO₃ of height `thickness` on top at an angle `ϕ`. The default for `ϕ` is chosen to be 45°, because that gives full conversion from s/p to p/s in the case of a λ/2-wave plate.
"""
TStructure(thickness, ϕ=45°) = LayeredStructure(
	superstrate = Si(),
	layers = [MoO₃(d = thickness, ϕ = ϕ)],
	substrate = Layer()
)

# ╔═╡ da9b9ff3-5dda-4515-922f-650051cf2327
"""
	RStructure(thickness, ϕ=45°)

Air as incident medium with α-MoO₃ of height `thickness` on top at of a gold mirror an angle `ϕ`. The default for `ϕ` is chosen to be 45°, because that gives full conversion from s/p to p/s in the case of a λ/2-wave plate.
"""
RStructure(thickness, ϕ=45°) = LayeredStructure(
	superstrate = Layer(),
	layers = [MoO₃(d = thickness, ϕ = ϕ)],
	substrate = Au()
)

# ╔═╡ ed602ad1-3bf4-463a-b14c-a3c871bd4c67
"""
	Rₚ(θ, n₁, n₂)

Reflection of s-polarised light at an interface with refractive indices `n₁` & `n₂` at an angle `θ`."""
function Rₚ(θ, n₁, n₂)
	num = n₁ * sqrt(1 - (n₁/n₂ * sin(θ))^2) - n₂ * cos(θ)
	denom = n₁ * sqrt(1 - (n₁/n₂ * sin(θ))^2) + n₂ * cos(θ)
	abs2(num / denom)
end

# ╔═╡ e34c1399-6a45-4597-898d-5956ec57731b
"""
Calculate transmittance form a array of StructureProperties. Returns an vectors with 4 arrays representing `Tₚₚ, Tₛₛ, Tₚₛ, Tₛₚ`. Where `f` is a prefactor needed in case the incidence and exit medium are not the same.
"""
function get_transmittance(Ps, f=1)
	T(x) = f * abs2(x)
	broadcast.(T, (unzip(transmission_coeffs.(Ps))))
end

# ╔═╡ fd3794cb-568a-4049-beda-9f6a8bf0e484
Tₚ_Air_Si = 1 - Rₚ(0, 1, n_Si)

# ╔═╡ 161c4991-035b-4923-abc9-de7383f4df1d
md"""
### Transmission
"""

# ╔═╡ 61358277-d585-4a0e-94a6-a192703b0e3d
begin #Simulation parameters
	d_wp_trans = (0.19:0.005:2.1)μm
	λ_wp_trans = (12.35:0.001:13.15)μm
	θ_wp_trans = 0°
	ζ_wp_trans = sin(θ_wp_trans)
end;

# ╔═╡ 0210c771-bb66-49ab-8562-e6c59fd59162
transmission_prefactor = 1/ n_Si

# ╔═╡ e224d40a-bc24-40ee-9728-4b1999b29abe
# ╠═╡ disabled = true
#=╠═╡
P_wp_trans = [calculate_structure_properties(ζ_wp_trans, λᵢ, TStructure(dᵢ))  
	for dᵢ ∈ d_wp_trans, λᵢ ∈ λ_wp_trans]
  ╠═╡ =#

# ╔═╡ c0af5a0e-5a13-4820-b123-131c78a8c5ba
#=╠═╡
T_wp_trans = Tₚ_Air_Si .* get_transmittance(P_wp_trans, transmission_prefactor);
  ╠═╡ =#

# ╔═╡ 30334532-5cd1-41bb-ad40-58e84a24b712
MK_T_wp_trans = Tₚ_Air_Si .* unzip([MacKay.calculate_transmission(λᵢ, θ_wp_trans, TStructure(dᵢ)) for dᵢ ∈ d_wp_trans, λᵢ ∈ λ_wp_trans])

# ╔═╡ 96b60f10-1c0c-4535-a310-a190d45b0e8a
#=╠═╡
let idx = 1
[heatmap(
		λ_wp_trans, d_wp_trans, T_wp_trans[idx],
		xlabel="Wavelength", ylabel="α-MoO₃ thickness", title="Transmittance",
		colorbar_title="Tₚₚ",
		xlims=(λ_wp_trans[1], λ_wp_trans[end]),
		ylims=(d_wp_trans[1], d_wp_trans[end]),
		clims=(0, 0.42), # x-phonon
), heatmap(
		λ_wp_trans, d_wp_trans, MK_T_wp_trans[idx],
		xlabel="Wavelength", ylabel="α-MoO₃ thickness", title="Transmittance",
		colorbar_title="Tₚₚ",
		xlims=(λ_wp_trans[1], λ_wp_trans[end]),
		ylims=(d_wp_trans[1], d_wp_trans[end]),
		clims=(0, 0.42), # x-phonon
)]
end
  ╠═╡ =#

# ╔═╡ d3bc111f-b8a6-4593-bae7-dd3d792d0d72
l = λ_wp_trans[1]

# ╔═╡ d0dd4af6-6688-46a2-96ed-582f6905c2fd
a = deg2rad.(0:0.1:30)

# ╔═╡ d8828061-bb7a-41d0-b0fa-ba2dcfd32b94
d = d_wp_trans[1]

# ╔═╡ 8bf9d1e3-140b-4d1b-9e30-75efee47d4b8
function calc_T(λ, α, S)
	cosa_out, n_inc, n_out, qₓ = MacKay.calculate_cos_αt(ustrip(m, λ), α, S)
	# transmission_prefactor = real(n_out * cos(a_out)) / real(n_inc * cos(α))

	P = calculate_structure_properties(real(qₓ), λ, S)

	r = reflection_coeffs(P)
	# t = transmission_coeffs(P)

	# T = transmission_prefactor .* abs2.(t)

	R = reflection(P)

	R
end

# ╔═╡ a7f5919a-a5e5-4a71-8bf2-6f53bf1d5736
function new_calc_T(λ, α, S)
	# t = MacKay.transmission_coeffs(ustrip(m, λ), α, S)
	r = MacKay.reflection_coeffs(ustrip(m,λ), α, S)

	# T = MacKay.calculate_transmission(λ, α, S)
	R = MacKay.calculate_reflection(λ, α, S)

	R
end

# ╔═╡ ff45ee50-6cfd-4878-8c04-6c0772dd0804
t = TStructure(d)

# ╔═╡ 4b06f029-7fa7-46c8-9476-1df95f2de6af
old = unzip(calc_T.(l, a, Ref(t)))

# ╔═╡ 11d4f800-efd7-4763-bbac-83c8cb6a1a22
new = unzip(new_calc_T.(l, a, Ref(t)))

# ╔═╡ ce1b4ed8-e2ed-421f-a861-05678c51ec15
let b = rad2deg.(a)
	plot(b, old[4], lw=4, label="old")
	plot!(b, new[4], lw=4, ls=:dash, label="new")
#plot!(b, R_an, lw=4, ls=:dash, label="analytical", c=:black)
	hline!([0, 1])
end

# ╔═╡ 49b5051d-93cb-42d2-aa2d-ba15ee4f6e71
md"""
### Analytical
"""

# ╔═╡ 3c310cba-9352-4846-9d22-1c28d3f51359
ϕ(λ, n, d) = 2π * n * d / λ

# ╔═╡ 784760db-4666-4397-b6c0-146042e15bb3
function r_stack(r₁₂, r₂₃, ϕ)
	num = r₁₂ + r₂₃ * exp(2im * ϕ)
	denom = 1 + r₁₂ * r₂₃ * exp(2im * ϕ)

	num / denom
end

# ╔═╡ 17f73a69-2a43-49e0-99bf-8d82bb1d53cf
n_Si

# ╔═╡ 122d3164-809e-4766-814d-0acde09ed14a
n_x = sqrt(ϵ_MoO₃(l)[1,1])

# ╔═╡ 8ddac059-f74f-4601-a742-07ad7fe683e1
r₁₂ = unzip(fresnel_coefficients_general.(n_Si, n_x, sin.(a)))

# ╔═╡ ea02d5d4-0a8d-43b6-aaac-ece7dc1370ef
sin_a2 = @. (n_Si / n_x) * sin(a)

# ╔═╡ a9a06cfc-9334-4e81-b8fa-2a697e509d66
r₂₃ = unzip(fresnel_coefficients_general.(n_x, 1.0, sin_a2))

# ╔═╡ 02c37bc2-4dc4-4548-8267-c54b75fbca88
R_an = real.(r_stack.(r₁₂[2], r₂₃[2], ϕ(l, n_x, d)))

# ╔═╡ f2e93d03-079e-43cb-aafd-895cd86b963b
function calculate_cos_αt(n1, n2, sin_αi)
	sin_αt = (n1 / n2) * sin_αi

	cos_αt = sqrt(1 - sin_αt^2)
	ret = real(n2 * cos_αt)
	imt = imag(n2 * cos_αt)
	if real(n2) > 0
		if imag(n2) > 0
			if imt < 0
			cos_αt *= -1
			end
		elseif isapprox(imag(n2),0,atol=1e-8)
			if isapprox(imt,0,atol=1e-14)
				if ret < 0
					cos_αt *= -1
				end
			elseif imt < 0
				cos_αt *= -1
			end
		end
	# Backwards traveling wave
	elseif real(n2) < 0 && imag(n2) < 0
		if imt < 0
			cos_αt *= -1
		end 
	end
	# This probably breaks for active media. The current approach would be to pray we choose the right branch.
	cos_αt
end

# ╔═╡ aa607576-84a5-46c2-a847-4f570a6285cd
uu = false ? 1 : 0

# ╔═╡ 7f4296c7-c6b8-4239-ba7d-78cbd16c7acf
"""
    fresnel_coefficients_general(n1, n2, θi)

Compute the Fresnel reflection coefficients for s-polarized (rs) and p-polarized (rp) light,
given complex refractive indices n1 and n2 (n1 = n1' + iκ1, n2 = n2' + iκ2)
and the angle of incidence θi (in radians).

Returns a tuple (rs, rp).
"""
function fresnel_coefficients_general(n1, n2, sin_θi)
    # Snell's law: sin(θt) = (n1 / n2) * sin(θi)
    sin_θt = (n1 / n2) * sin_θi
    
    # cos(θt) calculated using sqrt(1 - sin^2(θt)) for complex n1 and n2
    cos_θi = sqrt(1 - sin_θi^2)
	
    cos_θt = calculate_cos_αt(n1, n2, sin_θi)

    # Fresnel coefficients
    rs = (n1 * cos_θi - n2 * cos_θt) / (n1 * cos_θi + n2 * cos_θt)
    rp = (n2 * cos_θi - n1 * cos_θt) / (n2 * cos_θi + n1 * cos_θt)
    
    return (rs, rp)
end

# ╔═╡ 8825f80a-5023-4516-b652-9c801cf7171a
function test(n1, n2, sin_θi)
	# Snell's law: sin(θt) = (n1 / n2) * sin(θi)
    sin_θt = (n1 / n2) * sin_θi
    
    # cos(θt) calculated using sqrt(1 - sin^2(θt)) for complex n1 and n2
    cos_θi = sqrt(1 - sin_θi^2)
    cos_θt = sqrt(1 - sin_θt^2)
end

# ╔═╡ 28d68b9b-5e5f-41ec-95a8-4e832140dae9
md"""
### Reflection
"""

# ╔═╡ cd48d1c0-98f9-4a31-9036-22c99c4259f4
begin #Simulation parameters
	d_wp_ref = (0.19:0.005:2.1)μm
	λ_wp_ref = (12.4:0.005:15.1)μm
	θ_wp_ref = 0°
	ζ_wp_ref = sin(θ_wp_ref)
end;

# ╔═╡ 71a8107b-510f-41d4-bcdc-18f9f0906a14
# ╠═╡ disabled = true
#=╠═╡
P_wp_ref = [calculate_structure_properties(ζ_wp_ref, λᵢ, RStructure(dᵢ))  for dᵢ ∈ d_wp_ref, λᵢ ∈ λ_wp_ref];
  ╠═╡ =#

# ╔═╡ 6e613923-89d7-496b-ad2b-ac302239ed38
#=╠═╡
R_wp_ref = unzip(reflection.(P_wp_ref))
  ╠═╡ =#

# ╔═╡ 70528261-5b5e-452f-8f48-231bb3243fd8
# ╠═╡ disabled = true
#=╠═╡
MK_R_wp_ref = unzip([MacKay.calculate_reflection(λᵢ, θ_wp_ref, RStructure(dᵢ)) for dᵢ ∈ d_wp_ref, λᵢ ∈ λ_wp_ref]);
  ╠═╡ =#

# ╔═╡ d8d8fd8a-539e-4e29-85ef-77edfd323e0c
#=╠═╡
let idx = 3
	[heatmap(
		λ_wp_ref, d_wp_ref, R_wp_ref[idx],
		xlabel="Wavelength", ylabel="α-MoO₃ thickness", title="Reflectance",
		colorbar_title="Rₚₚ",
		xlims=(λ_wp_ref[1], λ_wp_ref[end]),
		ylims=(d_wp_ref[1], d_wp_ref[end]),
		clims=(0,0.6)
	),
	heatmap(
		λ_wp_ref, d_wp_ref, MK_R_wp_ref[idx],
		xlabel="Wavelength", ylabel="α-MoO₃ thickness", title="Reflectance",
		colorbar_title="Rₚₚ",
		xlims=(λ_wp_ref[1], λ_wp_ref[end]),
		ylims=(d_wp_ref[1], d_wp_ref[end]),
		clims=(0,0.6)
	)
	]
end
  ╠═╡ =#

# ╔═╡ 7067ca81-8146-4edb-b211-6f14789ef158
md"""
## Chirality paper
"""

# ╔═╡ 39176d99-eea0-4c0d-a3f4-3eb73d35ec5b
md"""
### Setup
"""

# ╔═╡ 5f6c21d9-daed-4f24-ba36-e6092969a07a
function ϵ_x_MoO₃_bridgman(λ)
    # Convert λ in meter to frequency in cm⁻¹
    f = 1 / (λ * 1e2)
    # Parameters from paper
    ϵ∞_x = 3.5
    fₗₒ_x₁, fₗₒ_x₂ = 534.3, 962 # [cm⁻¹]
    fₜₒ_x₁, fₜₒ_x₂ = 506.7, 810 # [cm⁻¹]
    γ_x₁, γ_x₂ = 49.1, 8.1# [cm⁻¹]

    return (ϵ∞_x * TMM.lorentz_osc(f, fₗₒ_x₁, fₜₒ_x₁, γ_x₁)
            * TMM.lorentz_osc(f, fₗₒ_x₂, fₜₒ_x₂, γ_x₂))
end

# ╔═╡ cd5194c9-234f-4a6e-afa3-393fed7cda79
function ϵ_y_MoO₃_bridgman(λ)
    # Convert λ in meter to frequency in cm⁻¹
    f = 1 / (λ * 1e2)
    # Parameters from paper
    ϵ∞_y = 4.0#4.87
    fₗₒ_y = 867#888.5 # [cm⁻¹]
    fₜₒ_y = 545.6#516.3 # [cm⁻¹]
    γ_y = 9.5#18.07 # [cm⁻¹]

    return (ϵ∞_y * TMM.lorentz_osc(f, fₗₒ_y, fₜₒ_y, γ_y))
end

# ╔═╡ fb057938-4e0f-47a3-a1b8-ce8e1784feb5
@permittivity "MoO₃_bridgman" λ -> Diagonal([ϵ_x_MoO₃_bridgman(λ),
                                           ϵ_y_MoO₃_bridgman(λ), TMM.ϵ_z_MoO₃(λ)]);

# ╔═╡ 923bbdcd-4bc7-48f9-9136-112a5db01345
CD(R, L) = R - L

# ╔═╡ b690f3ff-53a0-43cf-8f23-30b509edd869
S_CD_Twist4 = LayeredStructure(
	superstrate = Air,
	layers = [
	    MoO₃_bridgman(d = 0.846μm, ϕ = 90° - 25°),
	    MoO₃_bridgman(d = 0.803μm, ϕ = 14°)
	],
	substrate = Au()
)

# ╔═╡ 7fd95735-cfdd-487f-bd0f-71aac0380652
S_CD_Twist4_D = LayeredStructure(
	superstrate = Air,
	layers = [
	    MoO₃_bridgman(d = 0.846μm, ϕ = 90° - 25° + 45°),
	    MoO₃_bridgman(d = 0.803μm, ϕ = 14° + 45°)
	],
	substrate = Au()
)

# ╔═╡ c4b5575c-a355-4f10-81e5-8affa701129d
md"""
### Reflection
"""

# ╔═╡ 3579ed71-b5e6-404c-b72f-2f5846029302
λ_CD_ref = (10:0.001:15)μm

# ╔═╡ 715c7535-f8ee-47ba-8d09-85614036d517
reflection_coeffs(calculate_structure_properties(0,12.3μm, S_CD_Twist4))

# ╔═╡ a5cb32cd-1778-43b5-9a09-faf1f26e75d9
rₚₚ, rₛₛ, rₚₛ, rₛₚ = MacKay.reflection_coeffs(12.3e-6, 0, S_CD_Twist4)

# ╔═╡ 7434031a-373f-415d-a1b3-6b74b74de407
MacKay.reflection_coeffs(12.3e-6, 0, S_CD_Twist4; basis=:circular)

# ╔═╡ 0ee2d44e-f398-4000-be65-4ce01ef6a95a
r_RR = -((rₛₛ + rₚₚ) - im * (rₛₚ - rₚₛ)) / 2

# ╔═╡ b40458b6-539d-4d4e-873c-8a03c75a7e3d
r_LL = -((rₛₛ + rₚₚ) + im * (rₛₚ - rₚₛ)) / 2

# ╔═╡ 7be54ba4-56b6-4f97-8e1d-36def714fca2
r_LR = ((rₛₛ - rₚₚ) - im * (rₛₚ + rₚₛ)) / 2

# ╔═╡ eb369f71-1694-48c2-8b51-4da2bac57652
function MK_CD_ref(λ, S)
	R_RR, R_LL, R_RL, R_LR = MacKay.calculate_reflection(λ, 0, S; basis=:circular)

	R_R = R_RR + R_RL
    R_L = R_LL + R_LR

	CD(R_R, R_L)
end

# ╔═╡ 4eff2c09-1e3d-441b-813a-9836e35b2936
MK_CD_CD_ref = MK_CD_ref.(λ_CD_ref, Ref(S_CD_Twist4))

# ╔═╡ d4724f00-924c-4d5c-8bfb-6ae7528761b7
md"""
### Emission
"""

# ╔═╡ f27c190e-6953-42e4-9219-4c7760be0907
λ_CD_em = (12.25:0.01:13.5)μm

# ╔═╡ 26e43d85-7444-4480-9561-7f27a8874c92
basis_change(M, T) = inv(T) * M * T

# ╔═╡ dddad11e-2d1d-4c33-b545-dae78935c282
const T_mat = 1/sqrt(2) * [
	1 -im
	1 im
]

# ╔═╡ bb13bd90-6eab-4e79-8878-8bea1b6ed155
function circ_coeffs(cₚₚ, cₛₛ, cₚₛ, cₛₚ)
    c_linear = [
        cₚₚ cₛₚ
        cₚₛ cₛₛ
    ]

    c_RR, c_RL, c_LR, c_LL = basis_change(c_linear, inv(T_mat))
    return c_RR, c_LL, c_RL, c_LR
end

# ╔═╡ d03052db-aba3-4338-8c6f-1134dcfac6d0
function CD_ref(λ, S)
	P = calculate_structure_properties(0, λ, S)

	r = reflection_coeffs(P)
    R_RR, R_LL, R_RL, R_LR = abs2.(circ_coeffs(r...))

	R_R = R_RR + R_RL
    R_L = R_LL + R_LR

	CD(R_R, R_L)
end

# ╔═╡ 83915057-9aa4-45ca-a5c2-7f8cd6f545d0
CD_CD_ref = CD_ref.(λ_CD_ref, Ref(S_CD_Twist4))

# ╔═╡ 8d181466-a4bd-45e0-9d7d-70c349a4cfb8
plot(λ_CD_ref, CD_CD_ref);plot!(λ_CD_ref, MK_CD_CD_ref, ls=:dash)

# ╔═╡ 371b891b-2706-42f9-8a33-b6cd0a787e8d
function calculate_S(λ)
	P_H = calculate_structure_properties(0, λ, S_CD_Twist4)
	P_D = calculate_structure_properties(0, λ, S_CD_Twist4_D)

	Rₚₚ_H, Rₛₛ_H, Rₚₛ_H, Rₛₚ_H = reflection(P_H)
	Rₚₚ_D, Rₛₛ_D, Rₚₛ_D, Rₛₚ_D = reflection(P_D)
	
	r = reflection_coeffs(P_H)
	R_RR, R_LL, R_RL, R_LR = abs2.(circ_coeffs(r...))

	A_x = 1 - (Rₚₚ_H + Rₚₛ_H)
	A_y = 1 - (Rₛₛ_H + Rₛₚ_H)

	A_d = 1 - (Rₚₚ_D + Rₚₛ_D)
	A_a = 1 - (Rₛₛ_D + Rₛₚ_D)

	A_L = 1 - (R_RR + R_RL)
	A_R = 1 - (R_LL + R_LR)
	
	S₀ = (A_x + A_y) / 2	
	S₁ = (A_x - A_y) / 2
	S₂ = -(A_d - A_a) / 2
	S₃ = -(A_L - A_R) / 2

	S₀, S₁, S₂, S₃
end

# ╔═╡ dae11300-03d7-4dc1-adb0-751917ad67ee
function MK_calculate_S(λ)
	Rₚₚ_H, Rₛₛ_H, Rₚₛ_H, Rₛₚ_H = MacKay.calculate_reflection(λ, 0, S_CD_Twist4)
	Rₚₚ_D, Rₛₛ_D, Rₚₛ_D, Rₛₚ_D = MacKay.calculate_reflection(λ, 0, S_CD_Twist4_D)
	
	R_RR, R_LL, R_RL, R_LR = MacKay.calculate_reflection(λ, 0, S_CD_Twist4, basis=:circular)

	A_x = 1 - (Rₚₚ_H + Rₚₛ_H)
	A_y = 1 - (Rₛₛ_H + Rₛₚ_H)

	A_d = 1 - (Rₚₚ_D + Rₚₛ_D)
	A_a = 1 - (Rₛₛ_D + Rₛₚ_D)

	A_L = 1 - (R_RR + R_RL)
	A_R = 1 - (R_LL + R_LR)
	
	S₀ = (A_x + A_y) / 2	
	S₁ = (A_x - A_y) / 2
	S₂ = -(A_d - A_a) / 2
	S₃ = -(A_L - A_R) / 2

	S₀, S₁, S₂, S₃
end

# ╔═╡ f7389418-676c-4915-9ef3-c24dda24896d
S₀, S₁, S₂, S₃ = unzip(calculate_S.(λ_CD_em))

# ╔═╡ c1a11176-d965-4512-8abe-a67b5f3a6c58
MK_S₀, MK_S₁, MK_S₂, MK_S₃ = unzip(MK_calculate_S.(λ_CD_em))

# ╔═╡ 7f771e97-9f5d-4d52-8c45-d1b198e21404
plot(λ_CD_em, [S₀ S₁ S₂ S₃]); plot!(λ_CD_em, [MK_S₀ MK_S₁ MK_S₂ MK_S₃], ls=:dash)

# ╔═╡ 38cf1c83-9434-4496-86fa-34e9bc9463e9
md"""
## Passler, Paarman paper
"""

# ╔═╡ ec5b244a-89c3-4898-b107-8c23426c3a77
k_PP = 750:0.1:1050 #[cm⁻¹]

# ╔═╡ b37c9942-a0b5-4c08-9329-38689b5888de
λ_PP = 1e4 ./ k_PP * 1μm

# ╔═╡ 4fc30758-b28c-4695-b29b-09165933c2f1
function n_KRS5(λ)
	x = λ * 1e6
	sqrt(
		1 + 1.8293958./(1-0.0225./x.^2) + 
		1.6675593./(1-0.0625./x.^2) + 
		1.1210424./(1-0.1225./x.^2)+
		0.04513366./(1-0.2025./x.^2)+
		12.380234./(1-27089.737./x.^2)
	)
end

# ╔═╡ bed4c1fa-f340-4193-9705-2d02aaa4827d
@permittivity "KRS5" x -> n_KRS5(x)^2 * Diagonal(ones(3));

# ╔═╡ 10648535-b350-4c29-979f-de2fb865fba6
S_PP(d_air) = LayeredStructure(
	superstrate= KRS5(),
	layers = [Layer(d=d_air)],
	substrate=SiC()
)

# ╔═╡ f4797be4-8e92-481d-9516-e30dc4970bf9
airgap = [2, 3.5, 5.5, 7.5] .* 1e-6

# ╔═╡ 716e4640-7040-45c3-8b92-cec556304c07
function f(gap)
	R_PP, _, _, _ = unzip(MacKay.calculate_reflection.(λ_PP, 30°, Ref(S_PP(gap))))
	R_PP
end

# ╔═╡ 86c3623d-9e0d-4dac-afb0-1f56dadf831d
R_PP = f.(airgap)

# ╔═╡ 8f798a20-0659-4b33-a837-17dc0e839b41
plot(k_PP, R_PP)

# ╔═╡ Cell order:
# ╟─107c43c8-6955-46b9-9849-4d45ad308535
# ╠═686524a8-d6ed-4ee9-9f27-fb1502b613d1
# ╠═68f96782-42d1-470f-bf10-6eefcd77186c
# ╠═70e0c156-ce29-4927-94c7-64769c0261a0
# ╠═5e39aff5-7e98-4f68-aed5-a81ace47714f
# ╠═1f69912d-4737-46f9-878a-562d02828ebb
# ╠═92640359-3f64-4b9c-ac7c-97aa58cd6478
# ╠═984350af-8423-4238-8742-9e5a78293693
# ╠═214f6ace-39d3-4ad7-a05b-415fef8ca763
# ╠═69731a15-62c7-4a93-b183-ce94f2132a37
# ╟─5f2a7902-8cfa-45a0-8e36-172fd3542ba6
# ╟─94acbf39-b4d7-47df-b121-95292a75b3a5
# ╠═a0d2f37c-b411-4ee9-9913-e96204e0df64
# ╠═78a20169-1ae5-4175-b139-9e11adf8f6d9
# ╠═8193e9a3-1544-4bbd-a52a-7d4d7e47b769
# ╠═94464fa0-21b1-4e4b-bfcb-98158d53e0a0
# ╠═da9b9ff3-5dda-4515-922f-650051cf2327
# ╠═ed602ad1-3bf4-463a-b14c-a3c871bd4c67
# ╠═e34c1399-6a45-4597-898d-5956ec57731b
# ╠═fd3794cb-568a-4049-beda-9f6a8bf0e484
# ╟─161c4991-035b-4923-abc9-de7383f4df1d
# ╠═61358277-d585-4a0e-94a6-a192703b0e3d
# ╠═0210c771-bb66-49ab-8562-e6c59fd59162
# ╠═e224d40a-bc24-40ee-9728-4b1999b29abe
# ╠═c0af5a0e-5a13-4820-b123-131c78a8c5ba
# ╠═30334532-5cd1-41bb-ad40-58e84a24b712
# ╠═96b60f10-1c0c-4535-a310-a190d45b0e8a
# ╠═d3bc111f-b8a6-4593-bae7-dd3d792d0d72
# ╠═d0dd4af6-6688-46a2-96ed-582f6905c2fd
# ╠═d8828061-bb7a-41d0-b0fa-ba2dcfd32b94
# ╠═8bf9d1e3-140b-4d1b-9e30-75efee47d4b8
# ╠═a7f5919a-a5e5-4a71-8bf2-6f53bf1d5736
# ╠═ff45ee50-6cfd-4878-8c04-6c0772dd0804
# ╠═4b06f029-7fa7-46c8-9476-1df95f2de6af
# ╠═11d4f800-efd7-4763-bbac-83c8cb6a1a22
# ╠═ce1b4ed8-e2ed-421f-a861-05678c51ec15
# ╟─49b5051d-93cb-42d2-aa2d-ba15ee4f6e71
# ╠═3c310cba-9352-4846-9d22-1c28d3f51359
# ╠═02c37bc2-4dc4-4548-8267-c54b75fbca88
# ╠═784760db-4666-4397-b6c0-146042e15bb3
# ╠═17f73a69-2a43-49e0-99bf-8d82bb1d53cf
# ╠═122d3164-809e-4766-814d-0acde09ed14a
# ╠═8ddac059-f74f-4601-a742-07ad7fe683e1
# ╠═ea02d5d4-0a8d-43b6-aaac-ece7dc1370ef
# ╠═a9a06cfc-9334-4e81-b8fa-2a697e509d66
# ╠═f2e93d03-079e-43cb-aafd-895cd86b963b
# ╠═aa607576-84a5-46c2-a847-4f570a6285cd
# ╠═7f4296c7-c6b8-4239-ba7d-78cbd16c7acf
# ╠═8825f80a-5023-4516-b652-9c801cf7171a
# ╟─28d68b9b-5e5f-41ec-95a8-4e832140dae9
# ╠═cd48d1c0-98f9-4a31-9036-22c99c4259f4
# ╠═71a8107b-510f-41d4-bcdc-18f9f0906a14
# ╠═6e613923-89d7-496b-ad2b-ac302239ed38
# ╠═70528261-5b5e-452f-8f48-231bb3243fd8
# ╠═d8d8fd8a-539e-4e29-85ef-77edfd323e0c
# ╟─7067ca81-8146-4edb-b211-6f14789ef158
# ╟─39176d99-eea0-4c0d-a3f4-3eb73d35ec5b
# ╠═5f6c21d9-daed-4f24-ba36-e6092969a07a
# ╠═cd5194c9-234f-4a6e-afa3-393fed7cda79
# ╠═fb057938-4e0f-47a3-a1b8-ce8e1784feb5
# ╠═923bbdcd-4bc7-48f9-9136-112a5db01345
# ╠═b690f3ff-53a0-43cf-8f23-30b509edd869
# ╠═7fd95735-cfdd-487f-bd0f-71aac0380652
# ╟─c4b5575c-a355-4f10-81e5-8affa701129d
# ╠═3579ed71-b5e6-404c-b72f-2f5846029302
# ╠═715c7535-f8ee-47ba-8d09-85614036d517
# ╠═a5cb32cd-1778-43b5-9a09-faf1f26e75d9
# ╠═7434031a-373f-415d-a1b3-6b74b74de407
# ╠═0ee2d44e-f398-4000-be65-4ce01ef6a95a
# ╠═b40458b6-539d-4d4e-873c-8a03c75a7e3d
# ╠═7be54ba4-56b6-4f97-8e1d-36def714fca2
# ╠═d03052db-aba3-4338-8c6f-1134dcfac6d0
# ╠═eb369f71-1694-48c2-8b51-4da2bac57652
# ╠═83915057-9aa4-45ca-a5c2-7f8cd6f545d0
# ╠═4eff2c09-1e3d-441b-813a-9836e35b2936
# ╠═8d181466-a4bd-45e0-9d7d-70c349a4cfb8
# ╟─d4724f00-924c-4d5c-8bfb-6ae7528761b7
# ╠═f27c190e-6953-42e4-9219-4c7760be0907
# ╠═26e43d85-7444-4480-9561-7f27a8874c92
# ╠═dddad11e-2d1d-4c33-b545-dae78935c282
# ╠═bb13bd90-6eab-4e79-8878-8bea1b6ed155
# ╠═371b891b-2706-42f9-8a33-b6cd0a787e8d
# ╠═dae11300-03d7-4dc1-adb0-751917ad67ee
# ╠═f7389418-676c-4915-9ef3-c24dda24896d
# ╠═c1a11176-d965-4512-8abe-a67b5f3a6c58
# ╠═7f771e97-9f5d-4d52-8c45-d1b198e21404
# ╟─38cf1c83-9434-4496-86fa-34e9bc9463e9
# ╠═ec5b244a-89c3-4898-b107-8c23426c3a77
# ╠═b37c9942-a0b5-4c08-9329-38689b5888de
# ╠═4fc30758-b28c-4695-b29b-09165933c2f1
# ╠═bed4c1fa-f340-4193-9705-2d02aaa4827d
# ╠═10648535-b350-4c29-979f-de2fb865fba6
# ╠═f4797be4-8e92-481d-9516-e30dc4970bf9
# ╠═716e4640-7040-45c3-8b92-cec556304c07
# ╠═86c3623d-9e0d-4dac-afb0-1f56dadf831d
# ╠═8f798a20-0659-4b33-a837-17dc0e839b41
