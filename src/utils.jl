"""
    ppoints(n::Real, a=0.5)
    ppoints(x::AbstractVector, a=0.5)

Function to generate equi-spread probability values. The number of points `n` can be inferred from the length of the vector `x`.
"""
ppoints(n::Real; a = 1/2) = (n>0) ? ((1:n) .- a) / ( n .+ 1 .- 2a) : nothing
ppoints(x::AbstractVector; a = 1/2) = (length(x)>0) ? ((1:length(x)) .- a) / ( length(x) .+ 1 .- 2a) : nothing