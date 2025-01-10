# Code structure

```text
.
├── GeneralizedTransferMatrixMethod.jl
|   └──> Imports, exports, constants, references.
├── TMM.jl
|   |    Core matrices to construct the full transfer matrix,
|   └──> based on "Mackay, T. G. & Lakhtakia, A.  The Transfer-Matrix 
|        Method in Electromagnetics and Optics. vol. 1 (2020)."
├── OpticsFunctions.jl
|   |    Functions that calculate optical properties from 
|   └──> the transfer matrix.
├── Types.jl
|   └──> Composite types (Layer, LayeredStructure).
├── Permittivities.jl
|   └──> @permittivity macro and predefined premittivities.
└── HelperFunctions.jl
    └──> Helpers like Euler matrix, basis change, ...
```
