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
    expx=exp(x)
    if expx<0.99999
        return x/(expx-1)
    elseif expx < 1.00001
        return bernoulli_horner(x)
    else
        return x/(expx-1)
    end
end

##############################################################


const etruncmax = 100.0
const etruncmin = -100.0

function trexp(x)
  if(x<etruncmin)
    return 1.0/(exp(etruncmin)*(-x+etruncmin+1.0));
  elseif (x<etruncmax)
    return exp(x);
  else
      return exp(etruncmax)*(x-etruncmax+1.0);
  end
end


""" 
$(SIGNATURES)

Bernoulli function ``B(x)=\\frac{x}{e^x-1}`` for exponentially
fitted upwind, joint evaluation for positive and negative
argument

Usually, we need ``B(x), B(-x)`` togehter, 
and it is cheaper to calculate them together.

Returns two real numbers containing the result for argument
`x` and argument `-x`.

"""
function fbernoulli_pm(x)
    expx=exp(x)
    if abs(expx-1)>0.00001
        bp=x/(expx-1)
        bm=x/(1-1/expx)
        return bp,bm
    else
        y=bernoulli_horner(x)
        return y,x+y
    end
end
