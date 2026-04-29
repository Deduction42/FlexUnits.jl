
const UPREFERRED = Units{Dimensions{FixRat32}, AffineTransform{Float64}}[]

set_preferred_unit(u::Units; warn_failure=true) = set_preferred_unit!(UPREFERRED, u, warn_failure=warn_failure)

function set_preferred_unit!(unitset::AbstractVector{<:Units}, u::Units; warn_failure=false)
    du = dimensions(u)

    function abs_dims_match(x::Units)
        dx = dimension(x) 
        return (dx == du) || (dx == inv(du))
    end

    ii = findfirst(abs_dims_match, unitset)

    if isnothing(ii)
        push!(unitset, u)
        complexity_sort!(unitset)
    elseif dimensions(unitset[ii]) == dimensions(u)
        unitset[ii] = u  
    elseif warn_failure 
        @warn "Cannot set preferred unit $(u), inverse already exists"
    end

    return u 
end


function simplify(dref::D; unit_set=UPREFERRED) where D<:AbstractDimensions
    numerator = Pair{eltype(unit_set), Int64}[]
    denominator = Pair{eltype(unit_set), Int64}[]
    remainder = dref

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

    #return compound_unit(numerator, denominator, remainder)
    return (numerator, denominator, remainder) #Test this for now, and build compound_unit later
end

#Number of dimes dimensioon "dx" fits cleanly into dref
function clean_fit_div(remainder::AbstractDimensions, d::AbstractDimensions) 
    ii == 0
    while ii < 1000
        if (complexity(remainder) - complexity(remainder/d)) < complexity(d)
            return (ii, remainder/d) 
        else
            ii += 1
            d = d*d 
        end
    end
    error("Number of iterations exceeded 1000")
end

complexity_sort!(v::AbstractVector::Units) = sort!(v, by=complexity, rev=true)

function complexity(d::D) where D <: AbstractDimensions
    return sum(fn-> abs(getproperty(d, fn)), fieldnames(D))
end

complexity(d::StaticDims{D}) where D = complexity(D)
complexity(q::QuantUnion) = complexity(dimension(q))
complexity(u::Units) = complexity(dimension(u))