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
ζ = sin(α)

Properties  = calculate_structure_properties(ζ, λ, Structure)

@testset "TMM tests" begin
    @test reflection_coeffs(Properties) == (
        -0.10536220448702555 - 0.16854320092956956im,
        -0.141171513256365 - 0.21833196421280568im,
        0.00012224147659709362 - 9.236581349350117e-5im,
        -0.00012224147659687005 + 9.236581349360481e-5im)
    @test transmission_coeffs(Properties) == (
        -0.8310177856668711 + 0.5194977036114501im,
        -0.81085948958292 + 0.5242948352476914im,
        -0.0026032984617574036 - 0.004095179667185374im,
        -0.00260329846175737 - 0.004095179667185243im)
    @test reflection(Properties) == (0.03950800471395104,
                                     0.0675982427521139,
                                     2.3474422102954027e-8,
                                     2.3474422102918518e-8)
    @test transmission(ζ, Properties) == (0.9604919718116275,
                                          0.9324017337734648)
end
