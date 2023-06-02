#=
Form factor + edge velocity calculation
=#

################################################
"""
$(SIGNATURES)

Calculate node volume  and voronoi surface contributions for cell.
"""

function cellfactors!(T::Type{Edge1D}, ::Type{Cartesian1D}, coord, cellnodes, icell, nodefac, edgefac)
    en = local_celledgenodes(T)
    K = cellnodes[en[1, 1], icell]
    L = cellnodes[en[2, 1], icell]
    xK = coord[1, K]
    xL = coord[1, L]
    d = abs(xL[1] - xK[1])
    nodefac[1] = d / 2
    nodefac[2] = d / 2
    edgefac[1] = 1 / d
    nothing
end

function cellfactors!(T::Type{Edge1D}, ::Type{<:Polar1D}, coord, cellnodes, icell, nodefac, edgefac)
    en = local_celledgenodes(T)
    K = cellnodes[en[1, 1], icell]
    L = cellnodes[en[2, 1], icell]
    xK = coord[1, K]
    xL = coord[1, L]
    r0 = xK[1]
    r1 = xL[1]
    if r1 < r0
        r0 = xL[1]
        r1 = xK[1]
    end
    rhalf = 0.5 * (r1 + r0)
    # cpar[1]= π*(r1*r1-r0*r0);         # circular volume
    nodefac[1] = π * (rhalf * rhalf - r0 * r0)   # circular volume between midline and boundary
    nodefac[2] = π * (r1 * r1 - rhalf * rhalf)   # circular volume between midline and boundary
    edgefac[1] = 2.0 * π * rhalf / (r1 - r0)     # circular surface / width
    nothing
end

function cellfactors!(T::Type{Edge1D}, ::Type{<:Spherical1D}, coord, cellnodes, icell, nodefac, edgefac)
    en = local_celledgenodes(T)
    K = cellnodes[en[1, 1], icell]
    L = cellnodes[en[2, 1], icell]
    xK = coord[1, K]
    xL = coord[1, L]
    r0 = xK[1]
    r1 = xL[1]
    if r1 < r0
        r0 = xL[1]
        r1 = xK[1]
    end
    rhalf = 0.5 * (r1 + r0)
    nodefac[1] = π * (rhalf^3 - r0^3) * 4.0 / 3.0   # sphere volume between midline and boundary
    nodefac[2] = π * (r1^3 - rhalf^3) * 4.0 / 3.0   # sphere volume between midline and boundary
    edgefac[1] = 4.0 * π * rhalf^2 / (r1 - r0)      # circular surface / width
    nothing
end

function cellfactors!(T::Type{Triangle2D}, ::Type{Cartesian2D}, coord, cellnodes, icell, npar, epar)
    en = local_celledgenodes(T)
    @views n = cellnodes[:, icell]

    # Fill matrix of edge vectors
    V = @SMatrix[coord[1, n[en[1, 1]]]-coord[1, n[en[2, 1]]] coord[1, n[en[1, 2]]]-coord[1, n[en[2, 2]]] coord[1, n[en[1, 3]]]-coord[1, n[en[2, 3]]];
                 coord[2, n[en[1, 1]]]-coord[2, n[en[2, 1]]] coord[2, n[en[1, 2]]]-coord[2, n[en[2, 2]]] coord[2, n[en[1, 3]]]-coord[2, n[en[2, 3]]]]

    # Compute determinant 
    det = V[1, 3] * V[2, 2] - V[1, 2] * V[2, 3]
    vol = abs(0.5 * det)
    ivol = 1.0 / vol

    # squares of edge lengths
    dd = (V[1, 1] * V[1, 1] + V[2, 1] * V[2, 1],
          V[1, 2] * V[1, 2] + V[2, 2] * V[2, 2],
          V[1, 3] * V[1, 3] + V[2, 3] * V[2, 3])

    # contributions to \sigma_kl/h_kl
    epar[1] = (dd[2] + dd[3] - dd[1]) * 0.125 * ivol
    epar[2] = (dd[3] + dd[1] - dd[2]) * 0.125 * ivol
    epar[3] = (dd[1] + dd[2] - dd[3]) * 0.125 * ivol

    # contributions to \omega_k

    npar .= 0.0
    for i = 1:3
        npar[en[1, i]] += epar[i] * dd[i] * 0.25
        npar[en[2, i]] += epar[i] * dd[i] * 0.25
    end

    nothing
end

function cellfactors!(T::Type{Triangle2D}, ::Type{<:Cylindrical2D}, coord, cellnodes, icell, npar, epar)
    en = local_celledgenodes(T)

    @views n = cellnodes[:, icell]

    # Fill matrix of edge vectors
    V = @SMatrix[coord[1, n[en[1, 1]]]-coord[1, n[en[2, 1]]] coord[1, n[en[1, 2]]]-coord[1, n[en[2, 2]]] coord[1, n[en[1, 3]]]-coord[1, n[en[2, 3]]];
                 coord[2, n[en[1, 1]]]-coord[2, n[en[2, 1]]] coord[2, n[en[1, 2]]]-coord[2, n[en[2, 2]]] coord[2, n[en[1, 3]]]-coord[2, n[en[2, 3]]]]

    # Compute determinant 
    det = V[1, 3] * V[2, 2] - V[1, 2] * V[2, 3]
    area = abs(0.5 * det)

    # Integrate R over triangle (via quadrature rule)
    vol = 2.0 * π * area * (coord[1, n[1]] + coord[1, n[2]] + coord[1, n[3]]) / 3.0

    dd = (V[1, 1] * V[1, 1] + V[2, 1] * V[2, 1],
          V[1, 2] * V[1, 2] + V[2, 2] * V[2, 2],
          V[1, 3] * V[1, 3] + V[2, 3] * V[2, 3])

    emid = @SArray[0.5*(coord[1, n[en[1, 1]]] + coord[1, n[en[2, 1]]]) 0.5*(coord[1, n[en[1, 2]]] + coord[1, n[en[2, 2]]]) 0.5*(coord[1, n[en[1, 3]]] + coord[1, n[en[2, 3]]]);
                   0.5*(coord[2, n[en[1, 1]]] + coord[2, n[en[2, 1]]]) 0.5*(coord[2, n[en[1, 2]]] + coord[2, n[en[2, 2]]]) 0.5*(coord[2, n[en[1, 3]]] + coord[2, n[en[2, 3]]])]

    # use epar as temp storage
    @views tricircumcenter!(epar, coord[:, n[1]], coord[:, n[2]], coord[:, n[3]])
    cc = (epar[1], epar[2])

    r(p) = p[1]
    z(p) = p[2]

    for i = 1:3
        @views epar[i] = π * (r(cc) + r(emid[:, i])) * sqrt((r(cc) - r(emid[:, i]))^2 + (z(cc) - z(emid[:, i]))^2) / sqrt(dd[i])
    end

    function rintegrate(coord1::T1, coord2::T2, coord3::T3) where {T1, T2, T3}
        V11 = coord2[1] - coord1[1]
        V21 = coord2[2] - coord1[2]

        V12 = coord3[1] - coord1[1]
        V22 = coord3[2] - coord1[2]
        local det = V11 * V22 - V12 * V21
        a = abs(0.5 * det)
        2.0 * π * a * (r(coord1) + r(coord2) + r(coord3)) / 3.0
    end

    npar .= 0.0
    for i = 1:3
        @views coord1 = coord[:, n[en[1, i]]]
        @views coord2 = coord[:, n[en[2, i]]]

        @views npar[en[1, i]] += rintegrate(coord1, cc, emid[:, i])
        @views npar[en[2, i]] += rintegrate(coord2, cc, emid[:, i])
    end
    nothing
end

function cellfactors!(T::Type{Tetrahedron3D}, ::Type{Cartesian3D}, coord, cellnodes, icell, npar, epar)
    # Transferred from WIAS/pdelib, (c) J. Fuhrmann, H. Langmach, I. Schmelzer
    @views n = cellnodes[:, icell]
    en = local_celledgenodes(T)
    @views pp1 = en[1, :]
    @views pp2 = en[2, :]

    pi1 = (5, 6, 5, 1)
    pi2 = (6, 3, 1, 4)
    pi3 = (4, 2, 3, 2)

    po1 = (2, 1, 2, 6)
    po2 = (1, 4, 6, 3)
    po3 = (3, 5, 4, 5)

    # use epar as intermediate memory
    for i = 1:6
        p1 = cellnodes[pp1[i], icell]
        p2 = cellnodes[pp2[i], icell]
        dx = coord[1, p1] - coord[1, p2]
        dy = coord[2, p1] - coord[2, p2]
        dz = coord[3, p1] - coord[3, p2]
        epar[i] = dx * dx + dy * dy + dz * dz
    end

    # non-allocating tuple
    dd = (epar[1], epar[2], epar[3], epar[4], epar[5], epar[6])
    epar .= 0

    x1 = coord[1, n[2]] - coord[1, n[1]]
    y1 = coord[2, n[2]] - coord[2, n[1]]
    z1 = coord[3, n[2]] - coord[3, n[1]]

    x2 = coord[1, n[3]] - coord[1, n[1]]
    y2 = coord[2, n[3]] - coord[2, n[1]]
    z2 = coord[3, n[3]] - coord[3, n[1]]

    x3 = coord[1, n[4]] - coord[1, n[1]]
    y3 = coord[2, n[4]] - coord[2, n[1]]
    z3 = coord[3, n[4]] - coord[3, n[1]]

    det = (x1 * (y2 * z3 - y3 * z2) + x2 * (y3 * z1 - y1 * z3) + x3 * (y1 * z2 - y2 * z1))

    if det < 0
        det = -det
    end
    vol = det / 6
    vv = 96 * 6 * vol

    for i = 1:4
        npar[i] = 0.0
        i1 = pi1[i]
        i2 = pi2[i]
        i3 = pi3[i]

        h1 = dd[i1] * (dd[i2] + dd[i3] - dd[i1])
        h2 = dd[i2] * (dd[i3] + dd[i1] - dd[i2])
        h3 = dd[i3] * (dd[i1] + dd[i2] - dd[i3])

        df = h1 + h2 + h3
        vf = (h1 * dd[po1[i]] + h2 * dd[po2[i]] + h3 * dd[po3[i]] - 2 * dd[i1] * dd[i2] * dd[i3]) / (vv * df)
        epar[i1] += h1 * vf
        epar[i2] += h2 * vf
        epar[i3] += h3 * vf
    end

    for i = 1:6
        npar[pp1[i]] += epar[i]
        npar[pp2[i]] += epar[i]
        epar[i] = 6 * epar[i] / dd[i]
    end
    nothing
end

################################################
"""
$(SIGNATURES)

Calculate node volume  contributions for boundary face.
"""

function bfacefactors!(T::Type{Vertex0D}, ::Type{Cartesian1D}, coord, bfacenodes, ibface, nodefac, edgefac)
    nodefac[1] = 1.0
    nothing
end

function bfacefactors!(T::Type{Vertex0D}, ::Type{<:Polar1D}, coord, bfacenodes, ibface, nodefac, edgefac)
    inode = bfacenodes[1, ibface]
    r = coord[1, inode]
    nodefac[1] = 2 * π * r
    nothing
end

function bfacefactors!(T::Type{Vertex0D}, ::Type{<:Spherical1D}, coord, bfacenodes, ibface, nodefac, edgefac)
    inode = bfacenodes[1, ibface]
    r = coord[1, inode]
    nodefac[1] = 4 * π * r^2
    nothing
end

function bfacefactors!(T::Type{Edge1D}, ::Type{Cartesian2D}, coord, bfacenodes, ibface, nodefac, edgefac)
    en = local_celledgenodes(T)
    i1 = bfacenodes[en[1, 1], ibface]
    i2 = bfacenodes[en[2, 1], ibface]
    dx = coord[1, i1] - coord[1, i2]
    dy = coord[2, i1] - coord[2, i2]
    d = sqrt(dx * dx + dy * dy)
    nodefac[1] = d / 2
    nodefac[2] = d / 2
    edgefac[1] = 1 / d
    nothing
end

function bfacefactors!(T::Type{Edge1D}, ::Type{<:Cylindrical2D}, coord, bfacenodes, ibface, nodefac, edgefac)
    en = local_celledgenodes(T)
    i1 = bfacenodes[en[1, 1], ibface]
    i2 = bfacenodes[en[2, 1], ibface]
    r1 = coord[1, i1]
    r2 = coord[1, i2]
    z1 = coord[2, i1]
    z2 = coord[2, i2]
    dr = r1 - r2
    rmid = (r1 + r2) / 2
    dz = z1 - z2
    l = sqrt(dr * dr + dz * dz)
    nodefac[1] = π * (r1 + rmid) * l / 2
    nodefac[2] = π * (r2 + rmid) * l / 2
    nothing
end

function bfacefactors!(T::Type{Triangle2D}, ::Type{<:Cartesian3D}, coord, bfacenodes, ibface, npar, epar)
    # Transferred from WIAS/pdelib, (c) J. Fuhrmann, H. Langmach, I. Schmelzer

    en = local_celledgenodes(T)
    @views n = bfacenodes[:, ibface]

    epar .= 0
    for j = 1:3
        d = coord[j, n[en[1, 1]]] - coord[j, n[en[2, 1]]]
        epar[1] += d * d
        d = coord[j, n[en[1, 2]]] - coord[j, n[en[2, 2]]]
        epar[2] += d * d
        d = coord[j, n[en[1, 3]]] - coord[j, n[en[2, 3]]]
        epar[3] += d * d
    end

    dd = (epar[1], epar[2], epar[3])

    # Kanten-Flaechenanteile (ohne Abschneiden); epar als Hilfsfeld benutzt
    epar[1] = (dd[2] + dd[3] - dd[1]) * dd[1]
    epar[2] = (dd[3] + dd[1] - dd[2]) * dd[2]
    epar[3] = (dd[1] + dd[2] - dd[3]) * dd[3]
    vol = sqrt(epar[1] + epar[2] + epar[3]) * 0.25

    d = 1.0 / (8 * vol)

    # Knoten-Flaechenanteile (ohne Abschneiden)
    npar[1] = (epar[3] + epar[2]) * d * 0.25
    npar[2] = (epar[1] + epar[3]) * d * 0.25
    npar[3] = (epar[2] + epar[1]) * d * 0.25

    # Kantengewichte 
    epar[1] = epar[1] * d / dd[1]
    epar[2] = epar[2] * d / dd[2]
    epar[3] = epar[3] * d / dd[3]
    nothing
end

##################################################################

#
# TODO: this should be generalized for more quadrules
#
function integrate(coordl, coordr, hnormal, velofunc)
    wl = 1.0 / 6.0
    wm = 2.0 / 3.0
    wr = 1.0 / 6.0
    coordm = 0.5 * (coordl + coordr)
    (vxl, vyl) = velofunc(coordl[1], coordl[2])
    (vxm, vym) = velofunc(coordm[1], coordm[2])
    (vxr, vyr) = velofunc(coordr[1], coordr[2])
    return (wl * vxl + wm * vxm + wr * vxr) * hnormal[1] + (wl * vyl + wm * vym + wr * vyr) * hnormal[2]
end

"""
$(SIGNATURES)

Project velocity onto grid edges,
"""
function edgevelocities(grid, velofunc)
    @assert dim_space(grid) < 3

    cn = grid[CellNodes]
    ec = grid[EdgeCells]
    en = grid[EdgeNodes]
    coord = grid[Coordinates]

    velovec = zeros(Float64, num_edges(grid))
    if dim_space(grid) == 1
        for iedge = 1:num_edges(grid)
            K = en[1, iedge]
            L = en[2, iedge]
            elen = coord[1, L] - coord[1, K]
            vx, vy = velofunc((coord[1, K] + coord[1, L]) / 2, 0)
            velovec[iedge] = -elen * vx
        end
    else
        for iedge = 1:num_edges(grid)
            K = en[1, iedge]
            L = en[2, iedge]
            p1 = @MVector zeros(2)
            p2 = @MVector zeros(2)
            tricircumcenter!(p1,
                             coord[:, cn[1, ec[1, iedge]]],
                             coord[:, cn[2, ec[1, iedge]]],
                             coord[:, cn[3, ec[1, iedge]]])
            if ec[2, iedge] > 0
                tricircumcenter!(p2,
                                 coord[:, cn[1, ec[2, iedge]]],
                                 coord[:, cn[2, ec[2, iedge]]],
                                 coord[:, cn[3, ec[2, iedge]]])
            else
                p2 .= 0.5 * (coord[:, K] + coord[:, L])
            end
            hnormal = coord[:, K] - coord[:, L]
            velovec[iedge] = integrate(p1, p2, hnormal, velofunc)
        end
    end
    return velovec
end

"""
$(SIGNATURES)

Project velocity onto boundary face normals
"""
function bfacevelocities(grid, velofunc)
    @assert dim_space(grid) < 3
    bfacenodes = grid[BFaceNodes]
    coord = grid[Coordinates]
    bfacecells = grid[BFaceCells]
    bfacenormals = grid[BFaceNormals]
    bfr = grid[BFaceRegions]
    velovec = zeros(Float64, 2, num_bfaces(grid))
    if dim_space(grid) == 1
        for ibface = 1:num_bfaces(grid)
            vx, vy = velofunc(coord[1, bfacenodes[1, ibface]])
            velovec[ibface] = vx * bfacenormals[1, ibface]
        end
    else
        for ibface = 1:num_bfaces(grid)
            p1 = coord[:, bfacenodes[1, ibface]]
            p2 = coord[:, bfacenodes[2, ibface]]
            pm = 0.5 * (p1 + p2)
            velovec[1, ibface] = integrate(p1, pm, bfacenormals[:, ibface], velofunc)
            velovec[2, ibface] = integrate(pm, p2, bfacenormals[:, ibface], velofunc)
        end
    end
    return velovec
end

"""
$(SIGNATURES)

Node form factors per boundary face
"""
bfacenodefactors(sys) = sys.boundary_assembly_data.nodefactors
