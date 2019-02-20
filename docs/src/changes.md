# Changes


## V0.2, Feb 20, 2019

- Changed signature of all callback functions:
  This also allows to pass user defined arrays etc. to the callback functions.

  - Besides of `flux`, they now all have `node::TwoPointFluxFVM.Node`
    as a second argument.

  - `flux` has `edge::TwoPointFluxFVM.Edge` as a second argument

  - the `x` argument in `source!()` is omitted, the same data
     are now found in `node.coord`


- New method `edgelength(::TwoPointFluxFVM.Edge)`
  
## V0.1, Dec. 2018

- Initial release
