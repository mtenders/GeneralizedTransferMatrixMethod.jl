![header](./assets/banner.svg)

---

The ``4 \times 4`` transfer matrix method was initially introduced by Teitler
and Henvis [[1](#References)] and later popularized by Berreman
[[2](#References)]. This method provides a concise approach for computing the
reflection and transmission properties of anisotropic optical materials. By
employing a matrix-based formulation, it simplifies the calculation of these
parameters for any stack of planar layers and lends itself to numerical
implementation. The present package implements the generalized transfer matrix
formalism as presented in [[3-5](#References)]. The authors provide a Python
implementation ([pyGTM](https://github.com/pyMatJ/pyGTM)), as well as a Matlab
implementation ([First version (2019)](https://doi.org/10.5281/zenodo.601496),
[Updated version (2020)](https://zenodo.org/record/3648041)).

!!! note "Notation"

    The notation utilized in both the documentation and the source code closely
    adheres to the one introduced in [[3](#References)]. Any instances of
    deviation from this in the code are explicitly addressed within the
    corresponding docstring.

## Upcoming features

- Function to directly calculate the fields within the layers.
- Interface to easily read in tabulated permittivity data.

## References

1. [Teitler, S. & Henvis, B. W. Refraction in Stratified, Anisotropic
   Media*. J. Opt. Soc. Am., JOSA 60, 830–834
   (1970).](https://doi.org/10.1364/JOSA.60.000830) 
2. [Berreman, D. W. Optics in Stratified and Anisotropic Media: 4×4-Matrix
   Formulation. J. Opt. Soc. Am., JOSA 62, 502–510
   (1972).](https://doi.org/10.1364/JOSA.62.000502) 
3. [Passler, N. C. & Paarmann, A. Generalized 4 × 4 matrix formalism for light
   propagation in anisotropic stratified media: study of surface phonon
   polaritons in polar dielectric heterostructures. J. Opt. Soc. Am. B 34, 2128
   (2017)](http://doi.org/10.1364/JOSAB.34.002128). 
4. [Passler, N. C. & Paarmann, A. Generalized 4 × 4 matrix formalism for light
   propagation in anisotropic stratified media: study of surface phonon
   polaritons in polar dielectric heterostructures: erratum. J. Opt. Soc. Am. B
   36, 3246 (2019).](http://doi.org/10.1364/JOSAB.36.003246) 
5. [Passler, N. C., Jeannin, M. & Paarmann, A. Layer-resolved absorption of
   light in arbitrarily anisotropic heterostructures. Phys. Rev. B 101, 165425
   (2020).](http://doi.org/10.1103/PhysRevB.101.165425) 

