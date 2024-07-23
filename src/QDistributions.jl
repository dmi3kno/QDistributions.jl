"""
Placeholder for a short summary about QDistributions.
"""
module QDistributions

using Distributions

import Base: convert, rand, mod
import Distributions: partype, pdf, logpdf, logcdf, cdf, sampler, 
       @distr_support, @check_args, insupport, ContinuousUnivariateDistribution
import LogExpFunctions: log1mexp
import Random: AbstractRNG
import StatsBase: params
import Roots: fzero, Brent
import Statistics: quantile
import SpecialFunctions: zeta

include("FGLD.jl")
include("GLDcsw.jl")
include("utils.jl")

export FGLD, GLDcsw,
sqf, shape, location, scale, asymmetry, steepness, params, partype, 
cdf, ccdf, pdf, 
quantile, cquantile, qdf, dqf,
median, iqr, rskewness, rkurtosis, rand, ppoints

end # module
