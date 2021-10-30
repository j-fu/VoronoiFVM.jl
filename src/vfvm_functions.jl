##############################################################
"""
$(SIGNATURES)

Calculation of Bernoulli function via Horner scheme based on Taylor
coefficients around 0.
"""

function bernoulli_horner(x)
    y=x/47_900_160
    y=x*y
    y=x*(-1/1_209_600+y)
    y=x*y
    y=x*(1/30_240+y)
    y=x*y
    y=x*(-1/720+y)
    y=x*y
    y=x*(1/12+y)
    y=x*(-1/2+y)
    y=1+y
end


# Bernoulli thresholds optimized for Float64
const bernoulli_small_threshold=0.25
const bernoulli_large_threshold=40.0
##############################################################
"""
$(SIGNATURES)

Bernoulli function ``B(x)=\\frac{x}{e^x-1}`` for exponentially
fitted upwinding.

The name `fbernoulli` has been chosen to avoid confusion
with Bernoulli from JuliaStats/Distributions.jl

Returns a real number containing the result.
"""
function fbernoulli(x)
    if x<-bernoulli_large_threshold
        -x
    elseif x>bernoulli_large_threshold
	zero(x)
    else
        expx=exp(x)
        expxm1=expx-1.0
        if abs(expxm1)>bernoulli_small_threshold
            x/expxm1
        else
            bernoulli_horner(x)
        end
    end
end

##############################################################


""" 
$(SIGNATURES)

Bernoulli function ``B(x)=\\frac{x}{e^x-1}`` for exponentially
fitted upwind, joint evaluation for positive and negative
argument

Usually, we need ``B(x), B(-x)`` togehter, 
and it is cheaper to calculate them together.

Returns two real numbers containing the result for argument
`x` and argument `-x`.

The error in comparison with the evaluation of the original expression
with BigFloat is less than 1.0e-15
"""
function fbernoulli_pm(x)
    if x< -bernoulli_large_threshold
        return -x, zero(x)
    elseif x>bernoulli_large_threshold
	return zero(x),x
    else
        expx=exp(x)
        expxm1=expx-1.0
        if abs(expxm1)>bernoulli_small_threshold
            bp=x/expxm1
            bm=x/(1.0-1.0/expx)
            return bp,bm
        else
            y=bernoulli_horner(x)
            return y,x+y
        end
    end
end
