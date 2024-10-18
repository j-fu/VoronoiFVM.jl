# [`ExtendableFEMBase`](https://github.com/chmerdon/ExtendableFEMBase.jl/) Extension

The extension for [`ExtendableFEMBase`](https://github.com/chmerdon/ExtendableFEMBase.jl/) extends the funcitons below for solution types of [`ExtendableFEM`](https://github.com/chmerdon/ExtendableFEM.jl/).
The name of the extension originates from the solution type `FEVectorBlock` that is definded in `ExtendableFEMBase`.

```@docs
edgevelocities(grid, vel::FEVectorBlock; kwargs...)
bfacevelocities(grid, vel::FEVectorBlock; kwargs...)
```