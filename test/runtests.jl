using Unitful
using GeneralizedTransferMatrixMethod
using Test

using LinearAlgebra

nₒ(λ) = sqrt(1. +
    0.28604141 +
    1.07044083 / (1. - 1.00585997e-2 / λ^2) +
    1.10202242 / (1. - 100. / λ^2))
nₑ(λ) = sqrt(1. +
    0.28851804 +
    1.09509924 / (1. - 1.02101864e-2 / λ^2) +
    1.15662475 / (1. - 100. / λ^2))
perm_test(λ) = Diagonal([nₒ(λ), nₑ(λ), nₒ(λ)])

@permittivity "TestLayer" perm_test

@testset "Permittivity tests" begin
    λ = 3.543
    @test perm_test(λ) == ϵ_TestLayer(λ)
    @test perm_test(λ) == ϵ_TestLayer(λ * 1u"m")
    @test perm_test(λ) == TestLayer().ϵ(λ * 1e6u"μm")
end

Air = Layer()
Structure = LayeredStructure(
    superstrate = Air,
    layers = [TestLayer(d = 2.0u"μm", ϕ = 46.3u"°")],
    substrate = Air
)

λ = 643.0u"nm"
α = 74.3u"°"


@testset "TMM tests" begin
    @test all(reflection_coeffs(λ, α, Structure) .≈ (
        -0.10536220447767417 - 0.16854320092572742im,
        -0.14117151324375576 - 0.21833196420780057im,
        -0.00012224147660430226 + 9.236581346934042e-5im,
        0.00012224147660439114 - 9.236581346940777e-5im))
    @test all(transmission_coeffs(λ, α, Structure) .≈ (
        -0.8310177856836864 + 0.5194977035876823im,
        -0.810859489601413 + 0.524294835224868im,
        0.002603298461642966 + 0.004095179667275272im,
        0.0026032984616428958 + 0.004095179667275026im))
    @test all(calculate_reflection(λ, α, Structure) .≈ (0.03950800471395104,
                                                        0.0675982427521139,
                                                        2.3474422102954027e-8,
                                                        2.3474422102918518e-8))
    @test all(calculate_transmission(λ, α, Structure) .≈ (0.960468424155492,
                                                          0.9323781861201347,
                                                          2.3547659387657428e-5,
                                                          2.3547659387655046e-5))
end
