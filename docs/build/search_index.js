var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#TwoPointFluxFVM-1",
    "page": "Home",
    "title": "TwoPointFluxFVM",
    "category": "section",
    "text": "Solver for coupled nonlinear partial differential equations based on the two point flux finite volume method on admissible grids.This Julia package merges the ideas behind pdelib/fvsys with the multiple dispatch paradigm of Julia. It allows to avoid the implementation of function derivatives. It instead uses the ForwardDiff and DiffResults  to evaluate user fucntions along with their jacobians.So far this is merely a design study to learn and evaluate Julia.   It is however aimed to be feasible at least for small projects.It requires Julia 1.0.Documentation created with Documenter.jl resides here"
},

{
    "location": "#Typical-usage-1",
    "page": "Home",
    "title": "Typical usage",
    "category": "section",
    "text": "\"\"\"\nStructure containing  userdata information\n\"\"\"\nmutable struct Physics\n    reaction::Function\n    flux::Function\n    storage::Function\n    source::Function\n    eps::Float64 \n    Physics()=new()\nend\n\n\"\"\"\nReaction term\n\"\"\"\nphysics.reaction=function(physics,node,f,u)\n    f[1]=u[1]^2\nend\n\n\"\"\"\nStorage term\n\"\"\"\nphysics.storage=function(physics,node,f,u)\n    f[1]=u[1]\nend\n\n\"\"\"\nFlux term\n\"\"\"\nphysics.flux=function(physics,edge,f,uk,ul)\n    f[1]=this.eps*(uk[1]^2-ul[1]^2)\nend \n\n\n\"\"\"\nSource term\n\"\"\"\nphysics.source=function(physics,node,f)\n    f[1]=1.0e-4*node.coord[1]\nend \n\n"
},

{
    "location": "#[TwoPointFluxFVM.Grid](@ref)-1",
    "page": "Home",
    "title": "TwoPointFluxFVM.Grid",
    "category": "section",
    "text": "This is a simplex grid structure."
},

{
    "location": "#[TwoPointFluxFVM.Node](@ref)-1",
    "page": "Home",
    "title": "TwoPointFluxFVM.Node",
    "category": "section",
    "text": "This represents a node  in TwoPointFluxFVM.Grid."
},

{
    "location": "#[TwoPointFluxFVM.Edge](@ref)-1",
    "page": "Home",
    "title": "TwoPointFluxFVM.Edge",
    "category": "section",
    "text": "This represents an edge between two neigboring control volumes created from TwoPointFluxFVM.Grid.Currently, constructors are TwoPointFluxFVM.Grid(X::Array{Real,1}) for one-dimensional domains and TwoPointFluxFVM.Grid(X::Array{Float64,1},Y::Array{Float64,1}) for two-dimensional domains."
},

{
    "location": "#[TwoPointFluxFVM.System](@ref)-1",
    "page": "Home",
    "title": "TwoPointFluxFVM.System",
    "category": "section",
    "text": "From instances of  TwoPointFluxFVM.Graph and TwoPointFluxFVM.Physics, a TwoPointFluxFVM.System which contains all the necessary data for the solution of the nonlinear system described by them."
},

{
    "location": "install/#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "install/#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "So far, the package has not been registered with Julia.Steps:Create a directory for Julia packages which have not been registered, let us denote this as PKG_DIRcd to PKG_DIR and clone this repository:\nEither from WIAS rhodecode server     git clone  https://repos.wias-berlin.de/users/fuhrmann/projects/julia-packages/TwoPointFluxFVM Or from github:     git clone  https://github.com/j-fu/TwoPointFluxFVM.jl TwoPointFluxFVMAdd the following line to  the file .julia/config/startup.jl in your  home directory (create this file it does not exist). PKG_DIR must be the full path name of the directory.     push!(LOAD_PATH, \"PKG_DIR\")Now, import TwoPointFluxFVM should work in Julia scripts."
},

{
    "location": "changes/#",
    "page": "Changes",
    "title": "Changes",
    "category": "page",
    "text": ""
},

{
    "location": "changes/#Changes-1",
    "page": "Changes",
    "title": "Changes",
    "category": "section",
    "text": ""
},

{
    "location": "changes/#V0.3,-approaching-1",
    "page": "Changes",
    "title": "V0.3, approaching",
    "category": "section",
    "text": "Complete rewrite of assembly based on sparse matrix structure to store node_dof information.\nSolution is now a nnodes x nspecies sparse matrix\nThe wonderful array interface of Julia still provides slicing etc in oder to access particular species without need to write any bulk_solution stuff or whatever\nvalue() for debugging in physics functions"
},

{
    "location": "changes/#V0.2,-Feb-20,-2019-1",
    "page": "Changes",
    "title": "V0.2, Feb 20, 2019",
    "category": "section",
    "text": "Changed signature of all callback functions: This also allows to pass user defined arrays etc. to the callback functions. In particular, velocity vectors can be passed this way.\nBesides of flux!(), they now all have node::TwoPointFluxFVM.Node as a second argument.\nflux!() has edge::TwoPointFluxFVM.Edge as a second argument\nthe x argument in source!() is omitted, the same data  are now found in node.coordNew method edgelength(edge::TwoPointFluxFVM.Edge)"
},

{
    "location": "changes/#V0.1,-Dec.-2018-1",
    "page": "Changes",
    "title": "V0.1, Dec. 2018",
    "category": "section",
    "text": "Initial release"
},

{
    "location": "alldocs/#",
    "page": "API Documentation",
    "title": "API Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "alldocs/#TwoPointFluxFVM.SubGrid",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.SubGrid",
    "category": "type",
    "text": "struct SubGrid{Tc} <: AbstractGrid\n\nSubgrid of parent grid (mainly for visualization purposes). Intended to hold support of species which are not defined everywhere.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.SubGrid-Tuple{TwoPointFluxFVM.Grid,AbstractArray}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.SubGrid",
    "category": "method",
    "text": "function SubGrid(parent::Grid, \n                 subregions::AbstractArray; \n                 transform::Function=copytransform!,\n                 boundary=false)\n\nCreate subgrid of list of regions.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.add_boundary_species-Tuple{TwoPointFluxFVM.System,Integer,AbstractArray}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.add_boundary_species",
    "category": "method",
    "text": "function add_boundary_species(this::System, ispec::Integer, regions::AbstractArray)\n\nAdd species to a list of boundary regions. Species numbers for bulk and boundary species have to be distinct.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.add_species-Tuple{TwoPointFluxFVM.System,Integer,AbstractArray}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.add_species",
    "category": "method",
    "text": "function add_species(this::System,ispec::Integer, regions::AbstractArray)\n\nAdd species to a list of bulk regions. Species numbers for bulk and boundary species have to be distinct.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.cellmask!-Tuple{TwoPointFluxFVM.Grid,AbstractArray,AbstractArray,Integer}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.cellmask!",
    "category": "method",
    "text": "function cellmask!(grid::Grid,          \n                   maskmin::AbstractArray, # lower left corner\n                   maskmax::AbstractArray, # upper right corner\n                   ireg::Integer;          # new region number for elements under mask\n                   eps=1.0e-10)            # tolerance.\n\nEdit region numbers of grid cells via rectangular mask.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.dof-Union{Tuple{Tv}, Tuple{SysArray{Tv},Integer,Integer}} where Tv",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.dof",
    "category": "method",
    "text": "function dof(a::SysArray,ispec, inode)\n\nGet number of degree of freedom. Return 0 if species is not defined in node.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.fbernoulli-Tuple{Real}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.fbernoulli",
    "category": "method",
    "text": "function fbernoulli(x::Real)\n\nBernoulli function implementation for exponentially fitted finite volumes.\n\nThe name fbernoulli has been chosen to avoid confusion with Bernoulli from JuliaStats/Distributions.jl\n\nReturns a real number containing the result.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.fbernoulli_pm-Tuple{Real}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.fbernoulli_pm",
    "category": "method",
    "text": "function fbernoulli_pm(x::Real)\n\nBernoulli function implementation for exponentially fitted finite volumes, joint evaluation for positive and negative argument\n\nUsually, we need B(x), B(-x) togehter,  and it is cheaper to calculate them together.\n\nReturns two real numbers containing the result for argument x and argument -x.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.getdof-Tuple{TwoPointFluxFVM.SysArray,Integer}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.getdof",
    "category": "method",
    "text": "function getdof(a::SysArray,i::Integer)\n\nReturn  value for degree of freedom.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.integrate-Union{Tuple{Tv}, Tuple{System{Tv},Function,SysArray{Tv}}} where Tv",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.integrate",
    "category": "method",
    "text": "function integrate(this::System,F::Function,U)\n\nIntegrate solution vector over domain. Returns an Array{Int64,1} containing the integral for each species.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_bfaces-Tuple{TwoPointFluxFVM.Grid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_bfaces",
    "category": "method",
    "text": "num_bfaces(grid::Grid)\n\nNumber of boundary faces in grid.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_cells-Tuple{TwoPointFluxFVM.AbstractGrid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_cells",
    "category": "method",
    "text": "num_cells(grid)\n\nNumber of cells in grid\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_nodes-Tuple{TwoPointFluxFVM.AbstractGrid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_nodes",
    "category": "method",
    "text": "num_nodes(grid)\n\nNumber of nodes in grid\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_nodes-Tuple{TwoPointFluxFVM.SysArray}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_nodes",
    "category": "method",
    "text": "num_nodes(a::SysArray)\n\nNumber of nodes (size of second dimension) of solution array.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.setdof!-Tuple{TwoPointFluxFVM.SysArray,Any,Integer}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.setdof!",
    "category": "method",
    "text": "function setdof!(a::SysArray,v,i::Integer)\n\nSet value for degree of freedom.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.solve-Union{Tuple{Tv}, Tuple{System{Tv},SysArray{Tv}}} where Tv",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.solve",
    "category": "method",
    "text": "function solve(\n    this::System,            # Finite volume system\n    oldsol::Array{Tv,1};     # old time step solution resp. initial value\n    control=NewtonControl(), # Solver control information (optional)\n    tstep::Tv=Inf            # Time step size. Inf means  stationary solution. (optional)\n    )\n\nSolution method for instance of System.\n\nPerform solution of stationary system (if tstep==Inf) or implicit Euler time step system. \n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.unknowns-Union{Tuple{System{Tv}}, Tuple{Tv}} where Tv",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.unknowns",
    "category": "method",
    "text": "function unknowns(system)\n\nCreate a solution vector for system.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.value",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.value",
    "category": "function",
    "text": "value(x)\n\nExtract value from dual number. Use to debug physics callbacks. Re-exported from ForwardDiff.jl\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.AbstractGrid",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.AbstractGrid",
    "category": "type",
    "text": "abstract type AbstractGrid\n\nAbstract type for grid like datastructures.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.Edge",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.Edge",
    "category": "type",
    "text": "mutable struct Edge\n\nStructure holding local edge information.\n\nFields:\n\nindex::Int32\nregion::Int32\nnodeK::Int32\nnodeL::Int32\ncoordK::Array{Float64,1}\ncoordL::Array{Float64,1}\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.Grid",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.Grid",
    "category": "type",
    "text": "mutable struct Grid\n\nStructure holding grid data.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.Grid-Union{Tuple{Array{Tc,1}}, Tuple{Tc}} where Tc",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.Grid",
    "category": "method",
    "text": "Grid(X::Array{Tc,1})\n\nConstructor for 1D grid.\n\nConstruct 1D grid from an array of node cordinates. It creates two boundary regions with index 1 at the left end and index 2 at the right end.\n\nPrimal grid holding unknowns: marked by o, dual grid marking control volumes: marked by |.\n\n o-----o-----o-----o-----o-----o-----o-----o-----o\n |--|-----|-----|-----|-----|-----|-----|-----|--|\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.Grid-Union{Tuple{Tc}, Tuple{Array{Tc,1},Array{Tc,1}}} where Tc",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.Grid",
    "category": "method",
    "text": "Grid(X::Array{Tc,1},X::Array{Tc,1})\n\nConstructor for 2D grid from coordinate arrays.  Boundary region numbers count counterclockwise:\n\nlocation number\nsouth 1\neast 2\nnorth 3\nwest 4\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.NewtonControl",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.NewtonControl",
    "category": "type",
    "text": "Control parameter structure for Newton method.\n\nFields:\n\ntol_absolute::Float64 # Tolerance (in terms of norm of Newton update)\ntol_relative::Float64 # Tolerance (relative to the first residual)\ndamp_initial::Float64      # Initial damping parameter\ndamp_growth::Float64  # Damping parameter growth factor\nmax_iterations::Int32     # Maximum number of iterations\nmax_lureuse::Int32 # Maximum number of reuses of lu factorization\ntol_linear::Float64 # Tolerance of iterative linear solver\nverbose::Bool      # Verbosity\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.NewtonControl-Tuple{Any}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.NewtonControl",
    "category": "method",
    "text": "NewtonControl()\n\nDefault constructor\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.Node",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.Node",
    "category": "type",
    "text": "mutable struct Node\n\nStructure holding local node information. Fields:\n\nindex::Int32\nregion::Int32\ncoord::Array{Real,1}\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.SubgridSysArrayView",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.SubgridSysArrayView",
    "category": "type",
    "text": "struct SubgridSysArrayView{Tv} <: AbstractArray{Tv,2}\n\nStruct holding information for solution array view on subgrid\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.SysArray",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.SysArray",
    "category": "type",
    "text": "struct SysArray{Tv} <: AbstractArray{Tv,2}\n    node_dof::SparseMatrixCSC{Tv,Int16}\nend\n\nStruct holding solution information for system. Solution is stored in a sparse matrix structure.\n\nThis class plays well with the abstract array interface\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.System",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.System",
    "category": "type",
    "text": "mutable struct System{Tv}\n\nMain structure holding data for system solution.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.System-Tuple{TwoPointFluxFVM.Grid,Any,Integer}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.System",
    "category": "method",
    "text": "function  System(grid::Grid, physics::Any, maxspec::Integer)\n\nConstructor for System. physics provides some user data, maxspec is the maximum number of species.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#Base.copy-Union{Tuple{SysArray{Tv}}, Tuple{Tv}} where Tv",
    "page": "API Documentation",
    "title": "Base.copy",
    "category": "method",
    "text": "copy(this::SysArray)\n\nCreate a copy of solution array\n\n\n\n\n\n"
},

{
    "location": "alldocs/#Base.eltype-Tuple{TwoPointFluxFVM.AbstractGrid}",
    "page": "API Documentation",
    "title": "Base.eltype",
    "category": "method",
    "text": "eltype(grid)\n\nReturn element type of grid coordinates.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#Base.getindex-Tuple{TwoPointFluxFVM.SubgridSysArrayView,Integer,Integer}",
    "page": "API Documentation",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(aview::SubgridSysArrayView,ispec::Integer,inode::Integer)\n\nAccessor method for subgrid array view.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#Base.getindex-Tuple{TwoPointFluxFVM.SysArray,Integer,Integer}",
    "page": "API Documentation",
    "title": "Base.getindex",
    "category": "method",
    "text": " getindex!(a::SysArray, ispec, inode)\n\nAccessor for solution array.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#Base.setindex!-Tuple{TwoPointFluxFVM.SubgridSysArrayView,Any,Integer,Integer}",
    "page": "API Documentation",
    "title": "Base.setindex!",
    "category": "method",
    "text": "setindex!(aview::SubgridSysArrayView,v,ispec::Integer,inode::Integer)\n\nAccessor method for subgrid array view.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#Base.setindex!-Tuple{TwoPointFluxFVM.SysArray,Any,Integer,Integer}",
    "page": "API Documentation",
    "title": "Base.setindex!",
    "category": "method",
    "text": " setindex!(a::SysArray, v, ispec, inode)\n\nAccessor for solution array.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#Base.size-Tuple{TwoPointFluxFVM.SubgridSysArrayView}",
    "page": "API Documentation",
    "title": "Base.size",
    "category": "method",
    "text": "size(a::SubgridSysArrayView)\n\nReturn size of solution array view.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#Base.size-Tuple{TwoPointFluxFVM.SysArray}",
    "page": "API Documentation",
    "title": "Base.size",
    "category": "method",
    "text": "size(a::SysArray)\n\nReturn size of solution array.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#Base.view-Union{Tuple{Tv}, Tuple{SysArray{Tv},SubGrid}} where Tv",
    "page": "API Documentation",
    "title": "Base.view",
    "category": "method",
    "text": "view(a::SysArray{Tv},sg::SubGrid)\n\nCreate a view of the solution array on a subgrid.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.bfacefactors!-Tuple{TwoPointFluxFVM.Grid,Any,Any}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.bfacefactors!",
    "category": "method",
    "text": "bfacefactors!(grid::Grid,icell,nodefac)\n\nCalculate node volume  and voronoi surface contributions for boundary face.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.bfacenode-Tuple{TwoPointFluxFVM.Grid,Any,Any}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.bfacenode",
    "category": "method",
    "text": "bfacenode(grid::Grid,inode,ibface)\n\nIndex of boundary face node.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.celledgenode-Tuple{TwoPointFluxFVM.Grid,Any,Any,Any}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.celledgenode",
    "category": "method",
    "text": "celledgenode(grid::Grid,inode,iedge,icell)\n\nIndex of cell edge node.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.cellfactors!-Tuple{TwoPointFluxFVM.Grid,Any,Any,Any}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.cellfactors!",
    "category": "method",
    "text": "cellfactors!(grid::Grid,icell,nodefac,edgefac)\n\nCalculate node volume  and voronoi surface contributions for cell.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.cellnode-Tuple{TwoPointFluxFVM.AbstractGrid,Any,Any}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.cellnode",
    "category": "method",
    "text": "cellnode(grid,i,icell)\n\nReturn index of i-th local node in cell icell\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.dim_grid-Tuple{TwoPointFluxFVM.Grid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.dim_grid",
    "category": "method",
    "text": "dim_grid(grid)\n\nTopological dimension of grid\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.dim_space-Tuple{TwoPointFluxFVM.AbstractGrid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.dim_space",
    "category": "method",
    "text": "dim_space(grid)\n\nSpace dimension of grid\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.edgelength-Tuple{TwoPointFluxFVM.Edge}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.edgelength",
    "category": "method",
    "text": "function edgelength(edge::Edge)\n\nCalculate the length of an edge. \n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.geomspace",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.geomspace",
    "category": "function",
    "text": "function geomspace(a::Real, b::Real, ha::Real, hb::Real, tol=1.0e-10)\n\n(Try to) create a subdivision of interval (a,b) stored in the  returned array X such that \n\nX[1]==a, X[end]==b\n(X[2]-X[1])<=ha+tol*(b-a)\n(X[end]-X[end-1])<=hb+tol*(b-a)\nThere is a number q such that  X[i+1]-X[i] == q*(X[i]-X[i-1])\nX is the array with the minimal possible number of points with the above property\n\nCaveat: the algorithm behind this is  well tested but unproven.\n\nReturns an Array containing the points of the subdivision.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.inidirichlet!-Union{Tuple{Tv}, Tuple{System{Tv},SysArray{Tv}}} where Tv",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.inidirichlet!",
    "category": "method",
    "text": "function inidirichlet!(this::System,U)\n\nInitialize dirichlet boundary values for solution.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.is_boundary_species-Tuple{TwoPointFluxFVM.System,Integer}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.is_boundary_species",
    "category": "method",
    "text": "function is_boundary_species(this::System, ispec::Integer)\n\nCheck if species number corresponds to boundary species.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.is_bulk_species-Tuple{TwoPointFluxFVM.System,Integer}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.is_bulk_species",
    "category": "method",
    "text": "function is_bulk_species(this::System, ispec::Integer)\n\nCheck if species number corresponds bulk species.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.nodecoord-Tuple{TwoPointFluxFVM.AbstractGrid,Any}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.nodecoord",
    "category": "method",
    "text": "nodecoord(grid,inode)\n\nReturn view of coordinates of node inode.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_bfaceregions-Tuple{TwoPointFluxFVM.Grid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_bfaceregions",
    "category": "method",
    "text": "num_bfaceregions(grid::Grid)\n\nNumber of boundary face regions in grid.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_cellregions-Tuple{TwoPointFluxFVM.Grid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_cellregions",
    "category": "method",
    "text": "num_cellregions(grid::Grid)\n\nNumber of cell regions in grid.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_dof-Tuple{TwoPointFluxFVM.System}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_dof",
    "category": "method",
    "text": "num_dof(this::System)\n\nNumber of degrees of freedom for system.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_edges_per_cell-Tuple{TwoPointFluxFVM.Grid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_edges_per_cell",
    "category": "method",
    "text": "num_edges_per_cell(grid::Grid)\n\nNumber of edges per grid cell.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_nodes_per_bface-Tuple{TwoPointFluxFVM.Grid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_nodes_per_bface",
    "category": "method",
    "text": "num_nodes_per_bface(grid::Grid)\n\nNumber of nodes per boundary face\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_nodes_per_cell-Tuple{TwoPointFluxFVM.AbstractGrid}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_nodes_per_cell",
    "category": "method",
    "text": "num_nodes_per_cell(grid)\n\nReturn number of nodes per cell in grid.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_species-Tuple{TwoPointFluxFVM.SysArray}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_species",
    "category": "method",
    "text": "num_species(a::SysArray)\n\nNumber of species (size of first dimension) of solution array.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.num_species-Tuple{TwoPointFluxFVM.System}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.num_species",
    "category": "method",
    "text": "num_species(this::System)\n\nNumber of species in system\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.reg_bface-Tuple{TwoPointFluxFVM.Grid,Any}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.reg_bface",
    "category": "method",
    "text": "reg_bface(grid, ibface)\n\nBoundary region number for boundary face\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.reg_cell-Tuple{TwoPointFluxFVM.Grid,Any}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.reg_cell",
    "category": "method",
    "text": "reg_cell(grid,icell)\n\nBulk region number for cell\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.values-Tuple{TwoPointFluxFVM.SysArray}",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.values",
    "category": "method",
    "text": "values(a::SysArray)=a.node_dof\n\nArray of values in solution array.\n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.Dirichlet",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.Dirichlet",
    "category": "constant",
    "text": "Constant to be used as boundary condition factor  to mark Dirichlet boundary conditons.    \n\n\n\n\n\n"
},

{
    "location": "alldocs/#TwoPointFluxFVM.fbernoulli_eps",
    "page": "API Documentation",
    "title": "TwoPointFluxFVM.fbernoulli_eps",
    "category": "constant",
    "text": "Constant for switch between Taylor series and full implementation\n\n\n\n\n\n"
},

{
    "location": "alldocs/#API-Documentation-1",
    "page": "API Documentation",
    "title": "API Documentation",
    "category": "section",
    "text": "Modules = [TwoPointFluxFVM]\nOrder=[:type, :function, :macro, :constant]"
},

{
    "location": "allindex/#",
    "page": "API Index",
    "title": "API Index",
    "category": "page",
    "text": ""
},

{
    "location": "allindex/#API-Index-1",
    "page": "API Index",
    "title": "API Index",
    "category": "section",
    "text": "Pages   = [\"alldocs.md\"]\nModules = [TwoPointFluxFVM]\nOrder=[:type, :function, :macro, :constant]"
},

{
    "location": "examples/OneSpeciesNonlinearPoisson/#",
    "page": "1D Nonlinear Poisson equation with one species",
    "title": "1D Nonlinear Poisson equation with one species",
    "category": "page",
    "text": ""
},

{
    "location": "examples/OneSpeciesNonlinearPoisson/#D-Nonlinear-Poisson-equation-with-one-species-1",
    "page": "1D Nonlinear Poisson equation with one species",
    "title": "1D Nonlinear Poisson equation with one species",
    "category": "section",
    "text": "Solve the nonlinear Poisson equation-nabla varepsilon nabla u^2 + u^2 = 00001xin $\\Omega=(0,1)$ with boundary condition $u(0)=1$ and $u(1)=0.5$.using Markdown\nMarkdown.parse(\"\"\"\n```julia\n$(read(\"../../../examples/OneSpeciesNonlinearPoisson.jl\",String))\n```\n\"\"\")"
},

{
    "location": "examples/TwoSpeciesNonlinearPoisson/#",
    "page": "1D Nonlinear Poisson equation with two species",
    "title": "1D Nonlinear Poisson equation with two species",
    "category": "page",
    "text": ""
},

{
    "location": "examples/TwoSpeciesNonlinearPoisson/#D-Nonlinear-Poisson-equation-with-two-species-1",
    "page": "1D Nonlinear Poisson equation with two species",
    "title": "1D Nonlinear Poisson equation with two species",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse(\"\"\"\n```julia\n$(read(\"../../../examples/TwoSpeciesNonlinearPoisson.jl\",String))\n```\n\"\"\")"
},

{
    "location": "examples/IonicLiquid/#",
    "page": "1D Ionic Liquid",
    "title": "1D Ionic Liquid",
    "category": "page",
    "text": ""
},

{
    "location": "examples/IonicLiquid/#D-Ionic-Liquid-1",
    "page": "1D Ionic Liquid",
    "title": "1D Ionic Liquid",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse(\"\"\"\n```julia\n$(read(\"../../../examples/IonicLiquid.jl\",String))\n```\n\"\"\")"
},

{
    "location": "examples/NonlinearPoisson2D/#",
    "page": "2D Nonlinear Poisson equation",
    "title": "2D Nonlinear Poisson equation",
    "category": "page",
    "text": ""
},

{
    "location": "examples/NonlinearPoisson2D/#D-Nonlinear-Poisson-equation-1",
    "page": "2D Nonlinear Poisson equation",
    "title": "2D Nonlinear Poisson equation",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse(\"\"\"\n```julia\n$(read(\"../../../examples/NonlinearPoisson2D.jl\",String))\n```\n\"\"\")"
},

{
    "location": "examples/NonlinearPoisson2D_Reaction/#",
    "page": "2D Nonlinear Poisson equation with reaction",
    "title": "2D Nonlinear Poisson equation with reaction",
    "category": "page",
    "text": ""
},

{
    "location": "examples/NonlinearPoisson2D_Reaction/#D-Nonlinear-Poisson-equation-with-reaction-1",
    "page": "2D Nonlinear Poisson equation with reaction",
    "title": "2D Nonlinear Poisson equation with reaction",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse(\"\"\"\n```julia\n$(read(\"../../../examples/NonlinearPoisson2D_Reaction.jl\",String))\n```\n\"\"\")"
},

{
    "location": "examples/ThreeRegions1D/#",
    "page": "Differing species sets in regions, 1D",
    "title": "Differing species sets in regions, 1D",
    "category": "page",
    "text": ""
},

{
    "location": "examples/ThreeRegions1D/#Differing-species-sets-in-regions,-1D-1",
    "page": "Differing species sets in regions, 1D",
    "title": "Differing species sets in regions, 1D",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse(\"\"\"\n```julia\n$(read(\"../../../examples/ThreeRegions1D.jl\",String))\n```\n\"\"\")"
},

{
    "location": "examples/NonlinearPoisson2D_BoundaryReaction/#",
    "page": "2D Nonlinear Poisson equation with boundary reaction",
    "title": "2D Nonlinear Poisson equation with boundary reaction",
    "category": "page",
    "text": ""
},

{
    "location": "examples/NonlinearPoisson2D_BoundaryReaction/#D-Nonlinear-Poisson-equation-with-boundary-reaction-1",
    "page": "2D Nonlinear Poisson equation with boundary reaction",
    "title": "2D Nonlinear Poisson equation with boundary reaction",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse(\"\"\"\n```julia\n$(read(\"../../../examples/NonlinearPoisson2D_BoundaryReaction.jl\",String))\n```\n\"\"\")"
},

{
    "location": "examples/NonlinearPoisson1D_BoundarySpecies/#",
    "page": "1D two species system with boundary reaction and boundary species",
    "title": "1D two species system with boundary reaction and boundary species",
    "category": "page",
    "text": ""
},

{
    "location": "examples/NonlinearPoisson1D_BoundarySpecies/#D-two-species-system-with-boundary-reaction-and-boundary-species-1",
    "page": "1D two species system with boundary reaction and boundary species",
    "title": "1D two species system with boundary reaction and boundary species",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse(\"\"\"\n```julia\n$(read(\"../../../examples/NonlinearPoisson1D_BoundarySpecies.jl\",String))\n```\n\"\"\")"
},

{
    "location": "examples/NonlinearPoisson2D_BoundarySpecies/#",
    "page": "2D Nonlinear Poisson equation with boundary reaction and boundary species",
    "title": "2D Nonlinear Poisson equation with boundary reaction and boundary species",
    "category": "page",
    "text": ""
},

{
    "location": "examples/NonlinearPoisson2D_BoundarySpecies/#D-Nonlinear-Poisson-equation-with-boundary-reaction-and-boundary-species-1",
    "page": "2D Nonlinear Poisson equation with boundary reaction and boundary species",
    "title": "2D Nonlinear Poisson equation with boundary reaction and boundary species",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse(\"\"\"\n```julia\n$(read(\"../../../examples/NonlinearPoisson2D_BoundarySpecies.jl\",String))\n```\n\"\"\")"
},

]}
