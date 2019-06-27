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
function geomspace(a::Tv, b::Tv, ha::Tv, hb::Tv, tol=1.0e-10) where Tv
    
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

