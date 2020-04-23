# temporary replacement for corresponding XGrid stuff

abstract type AbstractElementGeometry end
abstract type AbstractElementGeometry0D <: AbstractElementGeometry end
abstract type AbstractElementGeometry1D <: AbstractElementGeometry end
abstract type AbstractElementGeometry2D <: AbstractElementGeometry end
abstract type AbstractElementGeometry3D <: AbstractElementGeometry end
abstract type AbstractElementGeometry4D <: AbstractElementGeometry end

abstract type Vertex0D <: AbstractElementGeometry0D end

abstract type Edge1D <: AbstractElementGeometry1D end

abstract type Polygon2D <: AbstractElementGeometry2D end
abstract type Triangle2D <: Polygon2D end
abstract type Quadrilateral2D <: Polygon2D end
abstract type Pentagon2D <: Polygon2D end
abstract type Hexagon2D <: Polygon2D end
abstract type Parallelogram2D <: Quadrilateral2D end

abstract type Circle2D <: AbstractElementGeometry3D end

abstract type Polyhedron3D <: AbstractElementGeometry end
abstract type Tetrahedron3D <: Polyhedron3D end
abstract type Hexahedron3D <: Polyhedron3D end
abstract type Parallelepiped3D <: Hexahedron3D end
abstract type Prism3D <: Polyhedron3D end
abstract type TrianglePrism3D <: Prism3D end
abstract type Sphere3D <: AbstractElementGeometry3D end

abstract type Polychoron4D <: AbstractElementGeometry4D end
abstract type HyperCube4D <: AbstractElementGeometry4D end





abstract type AbstractCoordinateSystem end
abstract type Cartesian1D   <: AbstractCoordinateSystem end    #x 
abstract type Cartesian2D   <: AbstractCoordinateSystem end    #x,y
abstract type Cartesian3D   <: AbstractCoordinateSystem end    #x,y,z
abstract type Cylindrical2D <: AbstractCoordinateSystem end  #r,z
abstract type Cylindrical3D <: AbstractCoordinateSystem end  #r,ϕ,z
abstract type Polar2D       <: AbstractCoordinateSystem end  #r,ϕ
abstract type Polar1D       <: AbstractCoordinateSystem end  #r (integral over ϕ)
abstract type Spherical3D   <: AbstractCoordinateSystem end  #r,ϕ,θ
abstract type Spherical1D   <: AbstractCoordinateSystem end  #r (integral over ϕ,\theta)

