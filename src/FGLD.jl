"""
    FGLD(α, β, δ, κ)

The *Flattened Skew-Logistic Distribution*, referred to as 
Flattened Generalized Logistic Distribution in Chakraborty & Sharma (2021) 
has location parameter `α`, scale `β > 0`, skewness parameter δ ∈ (0,1) and peak
flatness (steepness) parameter `κ > 0`. It is described by a quantile function

```math
S(u; α, β, δ, κ) = α + β\\left[(1-\\delta)\\ln(u)-\\delta\\ln(1-u)+\\kappa u\\right]
```
References: 

- Chakrabarty, T. K., & Sharma, D. (2018). The Quantile-Based Skew Logistic Distribution 
with Applications. In A. K. Chattopadhyay & G. Chattopadhyay (Eds.), Statistics 
and its Applications (pp. 51–73). Springer. 
https://doi.org/10.1007/978-981-13-1223-6_6

- Chakrabarty, T. K., & Sharma, D. (2021). A Generalization of the Quantile-Based 
Flattened Logistic Distribution. Annals of Data Science, 8(3), 603–627. 
https://doi.org/10.1007/s40745-021-00322-3

"""


struct FGLD{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    β::T
    δ::T
    κ::T
    FGLD{T}(α::T, β::T, δ::T, κ::T) where {T} = new{T}(α,β,δ,κ)
end

function FGLD(α::T, β::T, δ::T, κ::T; check_args::Bool=true) where {T <: Real}
    @check_args FGLD (β, β > zero(β)) (δ, zero(δ) <= δ <= one(δ)) (κ, κ >= zero(κ))
    return FGLD{T}(α,β,δ,κ)
end


FGLD(α::Real, β::Real, δ::Real, κ::Real; check_args::Bool=true) = FGLD(promote(α,β,δ,κ)...; check_args=check_args)
FGLD(α::Integer, β::Integer, δ::Integer, κ::Integer; check_args::Bool=true) = FGLD(float(α), float(β), float(δ), float(κ); check_args=check_args)
FGLD(α::Real=0.0) = FGLD(α, one(α), one(α)/2, one(α); check_args=false) # no need to check

@distr_support FGLD -Inf Inf

#### Conversions
function convert(::Type{FGLD{T}}, α::S, β::S, δ::S, κ::S) where {T<:Real, S<:Real}
    FGLD(T(α), T(β), T(δ), T(κ))
end
Base.convert(::Type{FGLD{T}}, d::FGLD) where {T<:Real} = FGLD{T}(T(d.α), T(d.β), T(d.δ), T(d.κ))
Base.convert(::Type{FGLD{T}}, d::FGLD{T}) where {T<:Real} = d

#### Parameters

location(d::FGLD) = d.α
scale(d::FGLD) = d.β
shape(d::FGLD) = (d.δ, d.κ)
asymmetry(d::FGLD) = d.δ
steepness(d::FGLD) = d.κ

params(d::FGLD) = (d.α, d.β, d.δ, d.κ)
@inline partype(d::FGLD{T}) where {T<:Real} = T


#### Statistics

mean(d::FGLD) = d.α + d.β*(2*d.δ-1+κ/2)
median(d::FGLD) = d.α+d.β*(log(2)*(1-2*d.δ)+d.κ/2)
iqr(d::FGLD) = d.β*(log(3)+d.κ/2)

# mode

function var(d::FGLD) 
    α, β, δ, κ = params(d)
    β^2*(1-4δ*(1-δ)+κ/2*(1+κ/6)+(π^2)/3*δ*(1-δ))
end
function skewness(d::FGLD) 
    α, β, δ, κ = params(d)
   (β^3)*(2*(-1+2δ*(3-2δ*(3-2δ)))+κ/4*(2δ-1)*(κ/3+3)-6δ*zeta(3)*(1-δ*(3-2δ)) )
end 
function kurtosis(d::FGLD)
    α, β, δ, κ = params(d)
    (β^4)*(9*(1-8δ*(1-δ*(3-2δ*(2-δ)))) + 9κ*(1/2-2δ*(1-δ)) + 2/9*(κ^2)*(5-11δ*(1-δ)) + (κ^3)/2*(1/3+κ/40)+ (δπ^4)/5*(4/3-δ*(13/3-3δ*(2-δ)))+ 2δπ^2*(1-δ*(5-4δ*(2-δ))) + κπ^2*(1+κ/6)*δ*(1-δ) )
end

#entropy

#### Functions

zval(d::FGLD, x::Real) = (x - d.α) / d.β
xval(d::FGLD, z::Real) = d.α + z * d.β

function quantile(d::FGLD, p::Real)
    α, β, δ, κ = params(d)
    α + β*((1-δ)*log(p)-δ*log(1-p)+κ*p)
end

function cquantile(d::FGLD, p::Real)
    α, β, δ, κ = params(d)
    α + β*((1-δ)*log(1-p)-δ*log(p)+κ*(1-p))
end

rskewness(d::FGLD) = (2*d.δ-1)*log(4/3)/(log(3)+d.κ/2) # Galton-Bowley skewness
rkurtosis(d::FGLD) = (log(21/5) + d.κ/2)/(log(3)+d.κ/2) # Moors' kurtosis

function lmoment(d::FGLD, r::Integer)
    α, β, δ, κ = params(d)
    if(r==1)
        α - β((1-2δ)+κ/2)
    elseif (r==2)
        β*(1/2+κ/6)
    else
        β*(2δ-1)^mod(r,2)/(r*(r-1))
    end
end

lvar(d::FGLD) = (3+d.κ)/(3*(2*d.α-2*(1-2*d.δ)+d.κ)) #L-coefficient of variation
lskewness(d::FGLD) = (2*d.δ-1)/(3+d.κ) #L-skewness ratio
lkurtosis(d::FGLD) = 1/(2*(3+d.κ)) #L-kurtosis ratio

qdf(d::FGLD, p::Real) = d.β*((1-d.δ)/p + d.δ/(1-p)+d.κ) # quantile density function
dqf(d::FGLD, p::Real) = p*(1-p)/(d.β*(1-d.δ+p*(2*d.δ-1)+d.κ*p*(1-p))) # density quantile function
logdqf(d::FGLD, p::Real) = log(p)+log(1-p)-log(d.β)-log(1-d.δ+p*(2*d.δ-1)+d.κ*p*(1-p))

function cdf(d::FGLD, x::Real)
    (insupport(d, x)) || DomainError("Data outside of the distribution support!")
    if (isinf(x)) return Float64(x>0) end # return zero or one if on the bound
    fzero(u -> quantile(d, u)-x, (0, 1), Brent()) 
end

ccdf(d::FGLD, x::Real) = 1-cdf(d, x) 
 
pdf(d::FGLD, x::Real) = dqf(d, cdf(d, x))
logpdf(d::FGLD, x::Real) = logdqf(d, cdf(d, x))


## Sampling
 
rand(rng::AbstractRNG, d::FGLD) = quantile(d, rand(rng))
sampler(rng::AbstractRNG, d::FGLD) = rand(rng::AbstractRNG, d::FGLD)
