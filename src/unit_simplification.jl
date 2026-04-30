using Revise 
using FlexUnits, .UnitRegistry
import FlexUnits: _pretty_unit_pwr, usymbol, uscale, dimtype

const UPREFERRED = Units{Dimensions{FixRat32}, AffineTransform{Float64}}[]

set_preferred_unit(u::Units; warn_failure=true) = set_preferred_unit!(UPREFERRED, u, warn_failure=warn_failure)

function set_preferred_unit!(unitset::AbstractVector{<:Units}, u::Units; warn_failure=true)
    du = dimension(assert_scalar(u))

    function abs_dims_match(x::Units)
        dx = dimension(x) 
        return (dx == du) || (dx == inv(du))
    end

    ii = findfirst(abs_dims_match, unitset)

    if isnothing(ii)
        push!(unitset, u)
        complexity_sort!(unitset)
    elseif dimension(unitset[ii]) == dimension(u)
        unitset[ii] = u  
    elseif warn_failure 
        @warn "Cannot set preferred unit $(u), inverse already exists"
    end

    return u 
end

function simplify(q::Quantity)
    u = convert(Units{dimtype(q), AffineTransform{Float64}}, simplify(dimension(q)))
    return uconvert(u, q)
end

simplify(dref::StaticDims{d}; unit_set=UPREFERRED) where d = simplify(d)

function simplify(d::D; unit_set=UPREFERRED) where D<:AbstractDimensions
    numerator = Pair{eltype(unit_set), Int64}[]
    denominator = Pair{eltype(unit_set), Int64}[]
    remainder = d

    for u in unit_set
        (n, remainder) = clean_fit_div(remainder, dimension(u))
        if (n > 0) 
            push!(numerator, u=>n)
            continue
        end

        (n, remainder) = clean_fit_div(remainder, inv(dimension(u)))
        if (n > 0)
            push!(denominator, u=>n)
        end
    end

    return compound_unit(numerator, denominator, remainder) #Test this for now, and build compound_unit later
end

function compound_unit(numerator::Vector{<:Pair{<:Units,<:Real}}, denominator::Vector{<:Pair{<:Units,<:Real}}, remainder::D) where D<:AbstractDimensions

    #Collect the remainder into numerator and denominator
    for fn in fieldnames(D)
        p = getproperty(remainder, fn)
        if p > 0
            d = D(; (fn => p,)...)
            push!(numerator, Units(d, AffineTransform(), Symbol(string(d))) => 1)
        elseif p < 0
            d = D(; (fn => -p,)...)
            push!(denominator, Units(d, AffineTransform(), Symbol(string(d))) => 1)
        end
    end

    numstr = if isempty(numerator)
        "1"
    else
        numstr_cat = join((_pretty_unit_pwr(p) for p in numerator), ' ')
        length(numerator) == 1 ? numstr_cat : "("*numstr_cat*")"
    end

    denstr = if isempty(denominator)
        "" 
    else
        denstr_cat = join((_pretty_unit_pwr(p) for p in denominator), ' ')
        "/"*(length(denominator) == 1 ? denstr_cat : "("*denstr_cat*")")
    end

    scale = prod(p->uscale(p[1]^p[2]), numerator, init=1.0) / prod(p->uscale(p[1]^p[2]), denominator, init=1.0)
    dims  = prod(p->dimension(p[1]^p[2]), numerator, init=ud"") / prod(p->dimension(p[1]^p[2]), denominator, init=ud"")

    return Units(dims, AffineTransform(scale=scale), Symbol(numstr*denstr))
end


#Number of dimes dimensioon "dx" fits cleanly into dref
function clean_fit_div(remainder::AbstractDimensions, d::AbstractDimensions) 
    ii = 0
    while ii < 1000
        Δcomplexity = complexity(remainder) - complexity(remainder/d)

        #If the complexity reduction is less than the complexity of d, stop here
        if Δcomplexity < complexity(d) 
            return (ii, remainder) 
        else
            ii += 1
            remainder = remainder/d
        end
    end
    error("Number of iterations exceeded 1000")
end

complexity_sort!(v::AbstractVector{<:Units}) = sort!(v, by=complexity, rev=true)

function complexity(d::D) where D <: AbstractDimensions
    return sum(fn-> abs(getproperty(d, fn)), fieldnames(D))
end

complexity(d::StaticDims{D}) where D = complexity(D)
complexity(q::QuantUnion) = complexity(dimension(q))
complexity(u::Units) = complexity(dimension(u))

_pretty_unit_pwr(p::Pair{<:Units, <:Real}) = _pretty_unit_pwr(usymbol(p[1]), p[2])

#================================================================================================================
# Test Code
================================================================================================================#

for u in [u"Ω", u"V", u"W", u"J", u"Pa", u"N", u"C", u"(m/s)"]
    set_preferred_unit(u)
end

q = simplify(1u"kg/L"*9.81u"m/s^2"*20u"cm")
q = simplify(5u"A^2"*0.01u"ohm")
q = simplify(5u"W"*sqrt(2u"V"))