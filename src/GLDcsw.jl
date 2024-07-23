"""
    GLDcsw(μ, σ, χ, ξ)

The *Generalized Lambda Distribution* (in Chalabi, Scott, Würtz parameterization) with location parameter `μ`, scale `σ > 0` and asymmetry χ ∈ (-1,1) 
and steepness ξ ∈ (0,1) has standard quantile function (sqf)

```math
S(u; χ, ξ) = 
\\begin{cases} 
\\ln(u)-\\ln(1-u), \\quad χ=0, ξ = 1/2 \\\
\\ln(u)-\\frac{1}{2α}\\left[(1-u)^{2 α} -1 \\right], \\quad χ \\neq 0, ξ = 1/2 (1+χ) \\\
\\frac{1}{2β}(u^{2β}-1)-\\ln(1-u), \\quad χ \\neq 0, ξ = 1/2 (1-χ) \\\
\\frac{1}{α+β}(u^{α + β} - 1)- \\frac{1}{α-β}\\left[(1-u)^{α-β}-1\\right], \\quad otherwise
\\end{cases} 
```

where 

```math
α = \\frac{1}{2}\\frac{\\frac{1}{2}-ξ}{\\sqrt{ξ(1-ξ)}} \\\\
β = \\frac{1}{2}\\frac{χ}{\\sqrt{1-χ^2}} 
```

References

- Chalabi, Y., Würtz, D. (2012). Flexible Distribution Modeling with the Generalized Lambda Distribution. ETH working paper. 2012-05. 
Online at https://www.rmetrics.org/WhitePapers
"""

struct GLDcsw{T<:Real} <: ContinuousUnivariateDistribution
    μ::T
    σ::T
    χ::T
    ξ::T
    GLDcsw{T}(μ::T, σ::T, χ::T, ξ::T) where {T} = new{T}(μ,σ,χ,ξ)
end

function GLDcsw(μ::T, σ::T, χ::T, ξ::T; check_args::Bool=true) where {T <: Real}
    @check_args GLDcsw (σ, σ > zero(σ)) (χ, -one(χ) < χ < one(χ)) (ξ, zero(ξ) < ξ < one(ξ))
    return GLDcsw{T}(μ,σ,χ,ξ)
end

GLDcsw(μ::Real, σ::Real, χ::Real, ξ::Real; check_args::Bool=true) = GLDcsw(promote(μ,σ,χ,ξ)...; check_args=check_args)
GLDcsw(μ::Integer, σ::Integer, χ::Integer, ξ::Integer; check_args::Bool=true) = GLDcsw(float(μ), float(σ), float(χ), float(ξ); check_args=check_args)
GLDcsw(μ::Real=0.0) = GLDcsw(μ, one(μ), zero(μ), one(μ)/2; check_args=false)

par_alpha(d::GLDcsw) = 0.5*(0.5-d.ξ)/sqrt(d.ξ*(1-d.ξ))
par_beta(d::GLDcsw) = 0.5*d.χ/sqrt(1-d.χ^2) 

sqf_lb(d::GLDcsw) = (d.ξ < (1+d.χ)/2) ? -inv(par_alpha(d) + par_beta(d)) : -Inf 
sqf_ub(d::GLDcsw) = (d.ξ < (1-d.χ)/2) ? inv(par_alpha(d) - par_beta(d)) : Inf 

# basic quantile function
function sqf(d::GLDcsw, p::Real) 
    α = par_alpha(d)
    β = par_beta(d)

    # lower bound
    if (p==0)  return sqf_lb(d)  end

    # upper bound
    if (p==1) return sqf_ub(d)  end

    # between bounds
    if (d.χ==0 && d.ξ == 0.5) #a=
        log(p) - log1p(-p)
    elseif ((d.χ != 0) && ((1+d.χ)/2-d.ξ == 0))
        log(p) - inv(2α)*((1-p)^(2α)-1)
    elseif ((d.χ != 0) && ((1-d.χ)/2-d.ξ == 0))
        inv(2β)*(p^(2β)-1) - log1p(-p)
    else
        inv(α+β)*(p^(α+β)-1) - inv(α-β)*((1-p)^(α-β)-1)
    end
end

sqf_iqr(d::GLDcsw) = sqf(d, 0.75)-sqf(d, 0.25)
sqf_median(d::GLDcsw) = sqf(d, 0.5)

# centered (and scaled) basic quantile function
csqf(d::GLDcsw, p::Real) = (sqf(d, p)- sqf_median(d))/sqf_iqr(d)

# lower and upper bound of the distribution
GLDcsw_lb(d::GLDcsw) = d.μ + d.σ*(sqf_lb(d) - sqf_median(d))/sqf_iqr(d)
GLDcsw_ub(d::GLDcsw) = d.μ + d.σ*(sqf_ub(d) - sqf_median(d))/sqf_iqr(d)

@distr_support GLDcsw GLDcsw_lb(d::GLDcsw) GLDcsw_ub(d::GLDcsw)


#### Conversions
function convert(::Type{GLDcsw{T}}, μ::S, σ::S, χ::S, ξ::S) where {T <: Real, S <: Real}
    GLDcsw(T(μ), T(σ), T(χ), T(ξ))
end
Base.convert(::Type{GLDcsw{T}}, d::GLDcsw) where {T<:Real} = GLDcsw{T}(T(d.μ), T(d.σ), T(d.χ), T(d.ξ))
Base.convert(::Type{GLDcsw{T}}, d::GLDcsw{T}) where {T<:Real} = d

#### Parameters

location(d::GLDcsw) = d.μ
scale(d::GLDcsw) = d.σ
shape(d::GLDcsw) = (d.χ, d.ξ)
asymmetry(d::GLDcsw) = d.χ
steepness(d::GLDcsw) = d.ξ

params(d::GLDcsw) = (d.μ, d.σ, d.χ, d.ξ)
@inline partype(d::GLDcsw{T}) where {T<:Real} = T

#### Statistics

median(d::GLDcsw) = d.μ
iqr(d::GLDcsw) = d.σ

#### Functions

zval(d::GLDcsw, x::Real) = (x - d.μ) / d.σ
xval(d::GLDcsw, z::Real) = d.μ + d.σ * z

# quantile function
quantile(d::GLDcsw, p::Real) =  d.μ + d.σ * csqf(d, p)
cquantile(d::GLDcsw, p::Real) = d.μ + d.σ * csqf(d, 1-p)

rskewness(d::GLDcsw) = (sqf(d, 0.75) + sqf(d,0.25) - 2*sqf_median(d))/sqf_iqr(d)
rkurtosis(d::GLDcsw) = (sqf(d, 7/8) - sqf(d,5/8) + sqf(d, 3/8) - sqf(d, 1/8))/sqf_iqr(d)

## lmoment
## lvar, lskewness, lkurtosis


# basic qdf
function sff(d::GLDcsw, p::Real) 
    α = par_alpha(d)
    β = par_beta(d)

    # between bounds
    if (d.χ==0 && d.ξ == 0.5) #both a=0, b=0
        inv(p*(1-p))       
    elseif ((d.χ != 0) && ((1+d.χ)/2-d.ξ == 0))
        inv(p) + (1-p)^(2α-1)
    elseif ((d.χ != 0) && ((1-d.χ)/2-d.ξ == 0))
        p^(2β-1) + inv(1-p)
    else
        p^(α+β-1) + (1-p)^(α-β-1)
    end
end

qdf(d::GLDcsw, p::Real) = d.σ/sqf_iqr(d)*sff(d,p)
dqf(d::GLDcsw, p::Real) = sqf_iqr(d)/(d.σ*sff(d,p))
logdqf(d::GLDcsw, p::Real) = log(dqf(d, p))

function cdf(d::GLDcsw, x::Real)
    (insupport(d, x)) || DomainError("Data outside of the distribution support!")
    if (x == GLDcsw_lb(d)) return 0 end
    if (x == GLDcsw_ub(d)) return 1 end
    Roots.fzero(u -> quantile(d, u)-x, (0, 1), Roots.Brent())
end
 
ccdf(d::GLDcsw, x::Real) = (1-cdf(d, x))

pdf(d::GLDcsw, x::Real) = dqf(d, cdf(d, x))
logpdf(d::GLDcsw, x::Real) = logdqf(d, cdf(d, x))

 ##### Sampling

rand(rng::AbstractRNG, d::GLDcsw) = quantile(d, rand(rng))
sampler(rng::AbstractRNG, d::GLDcsw) = rand(rng::AbstractRNG, d::GLDcsw)