# Changes

## V0.3, approaching

- Complete rewrite of assembly based on sparse matrix structure
  to store node_dof information.
  - Solution is now a nnodes x nspecies sparse matrix
  - The wonderful array interface of Julia still provides slicing
    etc in oder to access particular species without need to write
    any bulk_solution stuff or whatever

## V0.2, Feb 20, 2019

- Changed signature of all callback functions:
  This also allows to pass user defined arrays etc. to the callback functions.
  In particular, velocity vectors can be passed this way.

  - Besides of `flux!()`, they now all have `node::TwoPointFluxFVM.Node`
    as a second argument.

  - `flux!()` has `edge::TwoPointFluxFVM.Edge` as a second argument

  - the `x` argument in `source!()` is omitted, the same data
     are now found in `node.coord`


- New method `edgelength(edge::TwoPointFluxFVM.Edge)`
  
## V0.1, Dec. 2018

- Initial release
