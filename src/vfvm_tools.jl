##############################################################
"""
$(SIGNATURES)

(Try to) create a subdivision of interval (a,b) stored in the 
returned array X such that 
  - `X[1]==a, X[end]==b`
  - `(X[2]-X[1])<=ha+tol*(b-a)`
  - `(X[end]-X[end-1])<=hb+tol*(b-a)`
  - There is a number q such that  `X[i+1]-X[i] == q*(X[i]-X[i-1])`
  - X is the array with the minimal possible number of points with the above property
  
Caveat: the algorithm behind this is  well tested but unproven.

Returns an Array containing the points of the subdivision.
"""
function geomspace(a::Tv, b::Tv, ha::Tv, hb::Tv; tol=1.0e-10) where Tv
    
    function _geomspace0(l,h0, hl, tol=1.0e-10)
        
        @assert (l>0.0)
        @assert (h0>0.0)
        @assert (hl>=h0)
        @assert((hl+h0)<l)
        
        #  We need to  adjust two things:
        
        # The sum of the geometric progression must
        # match the length of the interval, so lmismatch
        # should be zero:
        function lmismatch(q,k)
            return l - h0*(1-q^k)/(1-q)
        end
        
        # The claim from experimenral evidence (Wolfram) 
        # is that, if written as a polynomial,
        # it has two real zeros: one and the value searched for which
        # is slightly larger than one. All other zeros are on one circle.
        
        # The size of the last interval should be close to 
        # to hl, so hmismatch should be close to one and not larger than one
        function  hmismatch(q,k)
            return  h0*q^(k-1)/hl
        end

        # define initial number of intervals from
        # average of minmal and maximal h
        n=Int32(ceil((2.0*l/(h0+hl))))
        
        if n==1
            n=2
        end
        
        # define initial q such that hmismatch is one.
        q=(hl/h0)^(1.0/(n-1.0))
        
        # Iteration until both mismatches are satisfactory
        # Outer loop runs until hmismatch is less than 1
        hmiss=10.0 # some initial value >1 just to run the loop at least once
        if abs(q-1.0)<tol
            hmiss=1.0
        end
        while  hmiss>1.0
            # increase number of intervals until
            # lmismatch becomes less than zero
            while  lmismatch(q,n)>0.0
                n+=1
            end
            
            # find initial interval for q containing
            # value with zero lmismatch 
            ns=0
            nsmax=1000
            
            while lmismatch(q,n)<0.0 &&  ns<nsmax
                q*=0.99
                ns+=1
            end
            
            xl=q
            xr=q/0.99
            @assert ns<nsmax
            
            # bisection to define q with zero lmismatch
            ns=0
            xm=0.5*(xl+xr)
            while (xr-xl)>tol && ns<nsmax
                ns+=1
                mmm=lmismatch(xm,n)
                if mmm==0.0
                    break
                elseif   lmismatch(xl,n)*mmm<0.0
                    xr=xm
                else
                    xl=xm
                end
                xm=0.5*(xl+xr)
            end
            q=xm
            @assert ns<nsmax
            hmiss=hmismatch(q,n)
            if hmiss>1.0 
                n=n+1
            end
        end
        #  printf("%d %g %g %g\n",n,q,lmismatch(q,n),hmismatch(q,n))
        
        X = Array{Tv,1}(undef,n+1)
        X[1]=0
        h=h0
        for i=1:n
            X[i+1]=X[i]+h
            h*=q
        end
        X[n+1]=l
        return X
    end

    # Map things to [0,b-a]
    @assert (ha>0.0)
    @assert (hb>0.0)
    @assert (a<b)


    tol=tol*(b-a)
    if ha<=hb
        X=_geomspace0(b-a,ha,hb,tol)
        X.+=a
    else
        X=-reverse(_geomspace0(b-a,hb,ha,tol))
        X.+=b
    end

    @assert (X[2]-X[1])<=ha+tol
    @assert (X[end]-X[end-1])<=hb+tol
    
    return X
end

"""
$(SIGNATURES)

Glue together two vectors a and b resulting in a vector c. They last element 
of a shall be equal (up to tol) to the first element of b.
The result fulfills `length(c)=length(a)+length(b)-1`
"""
function glue(a::Vector{Tv}, b::Vector{Tv}; tol=1.0e-10) where Tv
    #assert(is_monotone(a));
    #assert(is_monotone(b));
    na=length(a)
    nb=length(b)
    
    d=b[1]-a[na-1]
    @assert(d>0)
    d=b[1]-a[na]
    @assert(d>-tol)
    @assert(d<tol)

    c=Vector{Tv}(undef,na+nb-1)
    ic=0
    for ia=1:na
        ic+=1
        c[ic]=a[ia]
    end
    for ib=2:nb
        ic+=1
        c[ic]=b[ib]
    end
    return c
end


##################################################################
"""
$(SIGNATURES)

Extract value from dual number. Use to debug physics callbacks.
Re-exported from ForwardDiff.jl
"""
const value=ForwardDiff.value


##################################################################
"""
$(SIGNATURES)

Find the circumcenter of a triangle.                 
									 
Created from C source of Jonathan R Shewchuk <jrs@cs.cmu.edu>

Modified to return absolute coordinates.
"""
function tricircumcenter!(circumcenter::Vector{Tv},a::Vector{Tv},b::Vector{Tv},c::Vector{Tv}) where Tv

    # Use coordinates relative to point `a' of the triangle.
    xba = b[1] - a[1]
    yba = b[2] - a[2]
    xca = c[1] - a[1]
    yca = c[2] - a[2]

    # Squares of lengths of the edges incident to `a'.
    balength = xba * xba + yba * yba
    calength = xca * xca + yca * yca
    
    # Calculate the denominator of the formulae.
    # if EXACT
    #    Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html	 
    #    to ensure a correctly signed (and reasonably accurate) result, 
    #    avoiding any possibility of division by zero.		    
    #  denominator = 0.5 / orient2d((double*) b, (double*) c, (double*) a)
    
    
    # Take your chances with floating-point roundoff
    denominator = 0.5 / (xba * yca - yba * xca)

    # Calculate offset (from `a') of circumcenter. 
    xcirca = (yca * balength - yba * calength) * denominator  
    ycirca = (xba * calength - xca * balength) * denominator  


# The result is returned both in terms of x-y coordinates and xi-eta	 
# coordinates, relative to the triangle's point `a' (that is, `a' is	 
# the origin of both coordinate systems).	 Hence, the x-y coordinates	 
# returned are NOT absolute; one must add the coordinates of `a' to	 
# find the absolute coordinates of the circumcircle.  However, this means	 
# that the result is frequently more accurate than would be possible if	 
# absolute coordinates were returned, due to limited floating-point	 
# precision.  In general, the circumradius can be computed much more	 
# accurately.								 
    

circumcenter[1] = xcirca+a[1]
    circumcenter[2] = ycirca+a[2]

    return circumcenter
end




##################################################################
"""
$(SIGNATURES)

Project velocity onto grid edges,
"""
function edgevelocities(grid,velofunc)
    #
    # TODO: this should be generalized for more quadrules
    #
    function integrate(coordl,coordr,hnormal,velofunc)
        wl=1.0/6.0
        wm=2.0/3.0
        wr=1.0/6.0
        coordm=0.5*(coordl+coordr)
        (vxl,vyl)=velofunc(coordl[1],coordl[2])
        (vxm,vym)=velofunc(coordm[1],coordm[2])
        (vxr,vyr)=velofunc(coordr[1],coordr[2])
        return (wl*vxl +wm*vxm+wr*vxr)*hnormal[1] + (wl*vyl +wm*vym+wr*vyr)*hnormal[2]
    end
    
    @assert num_edges(grid)>0
    @assert dim_space(grid)<3

    cn=grid.cellnodes
    ec=grid.edgecells
    en=grid.edgenodes

    velovec=zeros(Float64,num_edges(grid))
    if VoronoiFVM.dim_space(grid)==1
        for iedge=1:num_edges(grid)
            K=en[1,iedge]
            L=en[2,iedge]
            elen=grid.coord[1,L]-grid.coord[1,K]
            vx,vy=velofunc((grid.coord[1,K]+grid.coord[1,L])/2,0)
            velovec[iedge]=-elen*vx
        end
    else
        for iedge=1:num_edges(grid)
            K=en[1,iedge]
            L=en[2,iedge]
            p1=Vector{Float64}(undef,2)
            p2=Vector{Float64}(undef,2)
            tricircumcenter!(p1,
                             grid.coord[:,cn[1,ec[1,iedge]]],
                             grid.coord[:,cn[2,ec[1,iedge]]],
                             grid.coord[:,cn[3,ec[1,iedge]]])
            if ec[2,iedge]>0
                tricircumcenter!(p2,
                                 grid.coord[:,cn[1,ec[2,iedge]]],
                                 grid.coord[:,cn[2,ec[2,iedge]]],
                                 grid.coord[:,cn[3,ec[2,iedge]]])
            else
                p2.=0.5*(grid.coord[:,K]+grid.coord[:,L])
            end    
            hnormal=grid.coord[:,K]-grid.coord[:,L]
            velovec[iedge]=integrate(p1,p2,hnormal,velofunc)
        end
    end
    return velovec
end
