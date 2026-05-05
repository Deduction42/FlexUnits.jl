#=============================================================================================
Broadcasting
=============================================================================================#

Base.BroadcastStyle(::Type{<:QuantUnion{T}}) where T = Base.BroadcastStyle(T)
Base.broadcastable(u::AbstractUnitLike) = Ref(u)

Base.size(q::QuantUnion) = size(ustrip(q))
Base.length(q::QuantUnion) = length(ustrip(q))
Base.axes(q::QuantUnion) = axes(ustrip(q))
Base.ndims(q::QuantUnion) = ndims(ustrip(q))
Base.ndims(::Type{<:QuantUnion{T}}) where {T} = ndims(T)
Base.iterate(q::Q, maybe_state...) where Q<:QuantUnion = 
    let subiterate = iterate(ustrip(q), maybe_state...)
        subiterate === nothing && return nothing
        return quantity(subiterate[1], unit(q)), subiterate[2]
    end

Base.size(u::AbstractUnitLike) = ()
Base.length(u::AbstractUnitLike) = 1
Base.axes(u::AbstractUnitLike) = ()
Base.ndims(u::AbstractUnitLike) = 0
Base.ndims(::Type{<:AbstractUnitLike}) = 0
Base.iterate(u::AbstractUnitLike) = (u, nothing)

Base.broadcastable(q::Q) where Q<:QuantUnion = quantity(Base.broadcastable(ustrip_base(q)), dimension(q))
Base.getindex(q::Q) where Q<:QuantUnion = quantity(getindex(ustrip(q)), unit(q))
Base.getindex(q::Q, inds) where Q<:QuantUnion = quantity(getindex(ustrip(q), inds), unit(q))
Base.getindex(q::Q, inds::CartesianIndex{0}) where Q<:QuantUnion = quantity(getindex(ustrip(q), inds), unit(q))
Base.getindex(q::Q, ind::Integer, inds::Integer...) where Q<:QuantUnion = quantity(getindex(ustrip(q), ind, inds...), unit(q))
Base.firstindex(q::Q) where Q<:QuantUnion = firstindex(ustrip(q))
Base.lastindex(q::Q) where Q<:QuantUnion = lastindex(ustrip(q))

#=============================================================================================
Displaying output
=============================================================================================#
#A global setting that can toggle "pretty printing" vs "parsable printing" 
#This value to "false" makes output parsable by default
const PRETTY_DIM_OUTPUT = Ref(true)
const DISPLAY_SIMPLIFIED = Ref(false)
const PREFERRED_UNITS = Units{Dimensions{FixRat32}, AffineTransform{Float64}}[]
const PREFERRED_LOGSCALE = Ref(dB)

function pretty_print_units(x::Bool) 
    PRETTY_DIM_OUTPUT[] = x 
    return x 
end

#Dispatching show methods
Base.show(io::IO, q::LogQuantUnion) = qshow(io, q, pretty=PRETTY_DIM_OUTPUT[])
Base.string(q::LogQuantUnion) = qstring(q, pretty=PRETTY_DIM_OUTPUT[])
Base.show(io::IO, u::AbstractUnitLike) = ushow(io, u, pretty=PRETTY_DIM_OUTPUT[])
Base.string(u::AbstractUnitLike) = ustring(u, pretty=PRETTY_DIM_OUTPUT[])

#=============================================================================================
Displaying quantities
=============================================================================================#
function qstring(x::LogQuantUnion; pretty=PRETTY_DIM_OUTPUT[])
    tmp_io = IOBuffer()
    qshow(tmp_io, x, pretty=pretty)
    return String(take!(tmp_io))
end

#Generic fallback for showing any uncaught quantity type
qshow(io::IO, q::LogQuantUnion; pretty=PRETTY_DIM_OUTPUT[]) = Base.show_default(io, q)

function qshow(io::IO, q::LogQuantUnion{<:Any, <:AbstractUnits}; pretty=PRETTY_DIM_OUTPUT[])
    if usymbol(unit(q)) != DEFAULT_USYMBOL #Print the unit symbol if it exists (not default)
        if pretty
            return qshow_pretty(io, q)
        else
            return qshow_parsable(io, q)
        end 

    elseif q isa LogQuant
        return qshow(io, logubase(q), pretty=pretty)

    else
        return qshow(io, ubase(q), pretty=pretty) #Print SI version (emphesizes that units are treated as SI quantities)
    end
end

function qshow(io::IO, q::QuantUnion{<:Any, <:AbstractDimLike}; pretty=PRETTY_DIM_OUTPUT[], simplified=DISPLAY_SIMPLIFIED[])
    if simplified
        return qshow(io, simplify(q))
    elseif pretty
        return qshow_pretty(io, q)
    else
        return qshow_parsable(io, q)
    end
end

function qshow(io::IO, q::LogQuant{<:Any, <:AbstractDimLike}; pretty=PRETTY_DIM_OUTPUT[], simplified=DISPLAY_SIMPLIFIED[])
    if simplified 
        return qshow(io, simplify(q))
    else
        print(io, "log(")
        qshow(io, ubase(q), pretty=pretty)
        return print(io, ")")
    end
end

function qshow_parsable(io::IO, q::LogQuantUnion)
    print(io, "(", ustrip(q), ")")
    return ushow(io, unit(q), pretty=false)
end

function qshow_pretty(io::IO, q::LogQuantUnion)
    print(io, ustrip(q), " ")
    return ushow(io, unit(q), pretty=true)
end

#=============================================================================================
Displaying units and dimensions
=============================================================================================#
function ustring(x::AbstractUnitLike; pretty=PRETTY_DIM_OUTPUT[])
    tmp_io = IOBuffer()
    ushow(tmp_io, x, pretty=pretty)
    return String(take!(tmp_io))
end

#Fallback for showing any other unit type 
ushow(io::IO, u::AbstractUnitLike; pretty=PRETTY_DIM_OUTPUT[]) = Base.show_default(io, u)

#Showing units (shows symbol by default, but falls back to Base.show_default)
function ushow(io::IO, u::AbstractUnits; pretty=PRETTY_DIM_OUTPUT[])
    usymb = usymbol(u)
    if usymb == DEFAULT_USYMBOL
        return Base.show_default(io, u)
    else
        return Base.print(io, usymb)
    end
end

#Fallback for showing dimensions with generic values
function ushow(io::IO, d::D; pretty=PRETTY_DIM_OUTPUT[]) where D<: AbstractDimensions
    vals = map(fn->fn=>getproperty(d,fn), dimension_names(D))
    return show(io, vals)
end

#StaticDims are shown exactly like dynamic ones
function ushow(io::IO, u::StaticDims; pretty=PRETTY_DIM_OUTPUT[])
    return ushow(io, udynamic(u), pretty=pretty)
end

#Showing dimensions with numeric values
function ushow(io::IO, d::D; pretty=PRETTY_DIM_OUTPUT[]) where D<:AbstractDimensions{<:Real}
    dimnames = collect(static_fieldnames(D))
    dimunits = unit_symbols(D)

    if isunknown(d)
        print(io, "?/?")
        return nothing 
    end

    usep = ifelse(pretty, " ", "*")
    abs_dim_pwr(dim::Symbol) = _unit_pwr_string(dimunits[dim], abs(d[dim]), pretty=pretty)

    function print_dim_series(series_io::IO, dims)
        if length(dims) == 1
            print(series_io, abs_dim_pwr(dims[1]))
        elseif length(dims) > 1
            series_str = "(" * join(abs_dim_pwr.(dims), usep) * ")"
            print(series_io, series_str)
        end
    end

    num = filter(k-> !ismissing(d[k]) && (d[k]>0), dimnames)
    den = filter(k-> !ismissing(d[k]) && (d[k]<0), dimnames)

    if isempty(num) & isempty(den)
    elseif isempty(den)
        print_dim_series(io, num)
    elseif isempty(num)
        print(io, "1/")
        print_dim_series(io, den)
    else
        print_dim_series(io, num)
        print(io, "/")
        print_dim_series(io, den)
    end

    return nothing
end

#=============================================================================================
Unit Simplification
=============================================================================================#

#Main dispatch unit to select preferred units, just in case dimensions are incompatible with PREFERRED_UNITS
preferred_units(::Type{<:Any}) = PREFERRED_UNITS

#Main methods exposed to the user
simplify(q::LogQuantUnion) = simplify(q, preferred_units(dimvaltype(q)))

set_preferred_unit(u::Units; warn_failure=true) = set_preferred_unit!(preferred_units(dimvaltype(u)), u, warn_failure=warn_failure)

function set_preferred_logscale(ls::LogScale)
    PREFERRED_LOGSCALE[] = ls 
    return ls 
end

function display_simplified_units(mode::Bool) 
    DISPLAY_SIMPLIFIED[] = mode 
    return mode
end

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
    else 
        warn_failure && @warn "Replacing $(unitset[ii]) with an inverse unit $(u)"
        unitset[ii] = u
    end

    return u 
end

#Helper methods methods exposed to the user
function simplify(q::QuantUnion, unit_set::AbstractVector)
    u = simplify(dimension(q), unit_set)
    return uconvert(u, q)
end

function simplify(q::LogQuant, unit_set::AbstractVector)
    u = PREFERRED_LOGSCALE[](simplify(dimension(q), unit_set))
    return uconvert(u, q)
end

function simplify(dref::StaticDims{d}, unit_set::AbstractVector) where d 
    u = simplify(d, unit_set)
    return Units{StaticDims{d}}(dimension(u), tobase(u), usymbol(u))
end

function simplify(d::AbstractDimensions, unit_set::AbstractVector)
    numerator = Pair{eltype(unit_set), Int64}[]
    denominator = Pair{eltype(unit_set), Int64}[]
    remainder = d

    for u in unit_set
        (n, remainder) = _clean_fit_div(remainder, dimension(u))
        if (n > 0) 
            push!(numerator, u=>n)
            continue
        end

        (n, remainder) = _clean_fit_div(remainder, inv(dimension(u)))
        if (n > 0)
            push!(denominator, u=>n)
        end
    end

    return compound_unit(numerator, denominator, remainder)
end

#Potentially useful to the user, not exported by default
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

    #Numerator as a string
    numstr = if isempty(numerator)
        ""
    else
        numstr_cat = join((_pretty_unit_pwr(p) for p in numerator), ' ')
        length(numerator) == 1 ? numstr_cat : "("*numstr_cat*")"
    end

    #Denominator as a string
    denstr = if isempty(denominator)
        "" 
    else
        denstr_cat = join((_pretty_unit_pwr(p) for p in denominator), ' ')
        prefix = isempty(numerator) ? "1/" : "/"
        prefix*(length(denominator) == 1 ? denstr_cat : "("*denstr_cat*")")
    end

    #Calculate the scale and dimension of the total unit
    scale = prod(p->uscale(p[1]^p[2]), numerator, init=1.0) / prod(p->uscale(p[1]^p[2]), denominator, init=1.0)
    dims  = prod(p->dimension(p[1]^p[2]), numerator, init=D()) / prod(p->dimension(p[1]^p[2]), denominator, init=D())

    return Units(dims, AffineTransform(scale=scale), Symbol(numstr*denstr))
end

#Measure the complexity of a dimension
function complexity(d::D) where D <: AbstractDimensions
    return sum(fn-> abs(getproperty(d, fn)), fieldnames(D))
end

complexity(d::StaticDims{D}) where D = complexity(D)
complexity(q::QuantUnion) = complexity(dimension(q))
complexity(u::Units) = complexity(dimension(u))

#Create the "ith" dimensional unit
function dimensional_unit(::Type{D}, di::Integer) where D <: AbstractDimensions
    f(ii) = ifelse(ii==di, true, false)
    n = length(fieldnames(D))

    if di > n
        throw(ArgumentError("Index argument '$(di)' is greather than the number of dimensions in $(D)"))
    else
        return D( ntuple(f, Val(n))... )
    end
end

#Number of dimes dimensioon "dx" fits cleanly into dref
function _clean_fit_div(remainder::AbstractDimensions, d::AbstractDimensions) 
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

_pretty_unit_pwr(p::Pair{<:Units, <:Real}) = _pretty_unit_pwr(usymbol(p[1]), p[2])


#=============================================================================================
Helper methods for showing numeric dimension powers
=============================================================================================#
function _unit_pwr_string(u::Symbol, p::Real; pretty=PRETTY_DIM_OUTPUT[])
    return pretty ? _pretty_unit_pwr(u, p) : _parsable_unit_pwr(u, p)
end

function _parsable_unit_pwr(u::Symbol, p::Real)
    if iszero(p)
        return ""
    elseif isone(p)
        return string(u)
    elseif isinteger(p)
        return string(u)*"^"*string(Int(p))
    else
        return string(u)*"^("*string(Rational(p))*")"
    end
end

function _pretty_unit_pwr(u::Symbol, p::Real)
    SUPERSCRIPT_MAPPING = ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')
    INTCHARS = ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')

    if iszero(p)
        return ""
    elseif isone(p)
        return string(u)
    else
        s = isinteger(p) ? string(Int(p)) : string(Rational(p))
        chars = map(replace(s, "//" => "ᐟ")) do c
            if c ∈ INTCHARS
                SUPERSCRIPT_MAPPING[parse(Int, c)+1]
            elseif c == '-'
                '⁻'
            else
                c
            end
        end
        return string(u)*join(chars)
    end
end
