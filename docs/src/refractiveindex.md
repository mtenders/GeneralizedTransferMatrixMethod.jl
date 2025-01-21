# RefractiveIndex.jl integration

The [@permittivity](@ref @permittivity) macro accepts the types defined by [RefractiveIndex.jl](https://github.com/stillyslalom/RefractiveIndex.jl).

After loading both packages
```@example refractiveindex
using GeneralizedTransferMatrixMethod
using RefractiveIndex
```
we can define some materials from [refractiveindex.info](https://refractiveindex.info/)
```@example refractiveindex
RI_Si = RefractiveMaterial("https://refractiveindex.info/?shelf=main&book=Si&page=Aspnes")
```
```@example refractiveindex
RI_SiO₂_o = RefractiveMaterial("https://refractiveindex.info/?shelf=main&book=SiO2&page=Ghosh-o")
```
```@example refractiveindex
RI_SiO₂_e = RefractiveMaterial("https://refractiveindex.info/?shelf=main&book=SiO2&page=Ghosh-e")
```
Now, we can use those together with the [@permittivity](@ref @permittivity) macro. For an isotropic material we can just pass one `RefractiveMaterial`.
```@example refractiveindex
@permittivity "Si" RI_Si

ϵ_Si(0.5e-6)
```
If we pass an array of `RefractiveMaterial` of length two, the second entry is used as the third diagonal component of the permittivity tensor.
```@example refractiveindex
@permittivity "SiO₂" [RI_SiO₂_o, RI_SiO₂_e]

ϵ_SiO₂(0.5e-6)
```
In all other cases, we can always just pass an array of length three.
```@example refractiveindex
@permittivity "SiO₂" [RI_SiO₂_o, RI_SiO₂_e, RI_SiO₂_o]

ϵ_SiO₂(0.5e-6)
```
!!! danger "TabulatedK"
    Materials with a table for just the extinction coefficient (`RefractiveMaterial{TabulatedK}`) are notd supported.
