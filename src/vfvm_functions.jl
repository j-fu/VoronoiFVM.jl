##############################################################
"""
$(SIGNATURES)


Constant for switch between Taylor series and full implementation
"""
const fbernoulli_eps=1.0e-4

##############################################################
"""
$(TYPEDSIGNATURES)

Bernoulli function implementation for exponentially
fitted finite volumes.



The name fbernoulli has been chosen to avoid confusion
with Bernoulli from JuliaStats/Distributions.jl

Returns a real number containing the result.
"""
function fbernoulli(x::Real)
    if x<-fbernoulli_eps
        return x/(exp(x)-1)
    elseif x < fbernoulli_eps
        x2  = x*x;
        x4  = x2*x2;
        x6  = x4*x2;
        x8  = x6*x2;
        x10 = x8*x2;
        return 1.0 - 0.5*x +1.0/12.0 * x2 
        - 1.0/720.0 * x4 
        + 1.0/30240.0 * x6 
        - 1.0/1209600.0 * x8 
        + 1.0/47900160.0 * x10; 
    else
        return x/(exp(x)-1)
    end
end

##############################################################
""" 
$(TYPEDSIGNATURES)

Bernoulli function implementation for exponentially
fitted finite volumes, joint evaluation for positive and negative
argument

Usually, we need B(x), B(-x) togehter, 
and it is cheaper to calculate them together.


Returns two real numbers containing the result for argument
`x` and argument `-x`.

"""
function fbernoulli_pm(x::Real)
    if x<-fbernoulli_eps
        expx=exp(x)
        bp=x/(exp(x)-1)
        bm=expx*bp
        return bp,bm
    elseif x < fbernoulli_eps
        x2  = x*x;
        x4  = x2*x2;
        x6  = x4*x2;
        x8  = x6*x2;
        x10 = x8*x2;

        horder=1.0/12.0 * x2 
        - 1.0/720.0 * x4 
        + 1.0/30240.0 * x6 
        - 1.0/1209600.0 * x8 
        + 1.0/47900160.0 * x10

        bp=  1.0 - 0.5*x + horder
        bm = 1.0 + 0.5*x + horder
        return bp,bm
    else
        expx=exp(x)
        bp=x/(exp(x)-1)
        bm=expx*bp
        return bp,bm
    end
end
