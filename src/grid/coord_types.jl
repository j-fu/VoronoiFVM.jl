##########################################################
"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 1D coordinates as cartesian.
"""
struct Cartesian1D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 2D coordinates as cartesian.
"""
struct Cartesian2D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 3D coordinates as cartesian.
"""
struct Cartesian3D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 1D coordinates as radial coordinate  r
assuming circular symmetry around origin 0.
"""
struct Polar1D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 2D coordinates as radial coordinate r
and eight coordinate z assuming circular symmetry 
around axis r=0.
"""
struct Cylindrical2D end

"""
$(TYPEDEF)

Type for dispatching formfactor calculations.
Intepret 1D coordinates as radial coordinate  r
assuming spehrical symmetry around point r=0.
"""
struct Spherical1D end


abstract type Simplex0D end
abstract type Simplex1D end
abstract type Simplex2D end
abstract type Simplex3D end
