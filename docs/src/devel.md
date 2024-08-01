# Development hints
Here, a bit of development hints are given which mainly concern tests and documentation generation.

## Pluto notebooks
The pluto notebooks in this package are "triple use":
- As typical Pluto notebooks, they are self-contained in the senses that they contain their own Project and Manifest files. So users can just download and execute them.
- If they run with the environmnet variable `PLUTO_PROJECT` set to some julia environment, this environment will activated at the start of the notebook. In particular, they use `Revise.jl` so they can be run during development of VoronoiFVM.jl. See also  https://github.com/fonsp/Pluto.jl/issues/1788 .
- During CI tests, they are run as scripts. For this purpose they are wrapped into temporary modules, and @test macros can be used in the notebooks.





