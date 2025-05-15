#=============================================================================================
Broadcasting
=============================================================================================#

Base.BroadcastStyle(::Type{<:AbstractQuantity{T}}) where T = Base.BroadcastStyle(T)
Base.broadcastable(u::AbstractUnitLike) = Ref(u)

Base.size(q::AbstractQuantity) = size(ustrip(q))
Base.length(q::AbstractQuantity) = length(ustrip(q))
Base.axes(q::AbstractQuantity) = axes(ustrip(q))
Base.ndims(q::AbstractQuantity) = ndims(ustrip(q))
Base.ndims(::Type{<:AbstractQuantity{T}}) where {T} = ndims(T)
Base.iterate(q::Q, maybe_state...) where Q<:AbstractQuantity = 
    let subiterate = iterate(ustrip(q), maybe_state...)
        subiterate === nothing && return nothing
        return constructorof(Q)(subiterate[1], unit(q)), subiterate[2]
    end

Base.size(u::AbstractUnitLike) = ()
Base.length(u::AbstractUnitLike) = 1
Base.axes(u::AbstractUnitLike) = ()
Base.ndims(u::AbstractUnitLike) = 0
Base.ndims(::Type{<:AbstractUnitLike}) = 0
Base.iterate(u::AbstractUnitLike) = (u, nothing)

Base.broadcastable(q::Q) where Q<:AbstractQuantity = Q(Base.broadcastable(ustrip_base(q)), dimension(q))
Base.getindex(q::Q) where Q<:AbstractQuantity = Q(getindex(ustrip(q)), unit(q))
Base.getindex(q::Q, inds) where Q<:AbstractQuantity = quantity(getindex(ustrip(q), inds), unit(q))
Base.getindex(q::Q, inds::CartesianIndex{0}) where Q<:AbstractQuantity = quantity(getindex(ustrip(q), inds), unit(q))
Base.getindex(q::Q, ind::Integer, inds::Integer...) where Q<:AbstractQuantity = quantity(getindex(ustrip(q), ind, inds...), unit(q))

#=============================================================================================
Displaying output
=============================================================================================#
#A global setting that can toggle "pretty printing" vs "parsable printing" 
#This value to "false" makes output parsable by default
const PRETTY_DIM_OUTPUT = Ref(true)

function Base.string(x::Union{AbstractUnitLike, AbstractQuantity}; pretty=PRETTY_DIM_OUTPUT[])
    tmp_io = IOBuffer()
    show(tmp_io, x, pretty=pretty)
    return String(take!(tmp_io))
end

function pretty_print_units(x::Bool) 
    PRETTY_DIM_OUTPUT[] = x 
    return x 
end

function Base.show(io::IO, q::AbstractQuantity{<:Any, <:AbstractDimensions}; pretty=PRETTY_DIM_OUTPUT[])
    if pretty
        return _pretty_print_quantity(io, q)
    else
        return _parsable_print_quantity(io, q)
    end
end

function Base.show(io::IO, q::AbstractQuantity{<:Any, <:AbstractUnits}; pretty=PRETTY_DIM_OUTPUT[])
    if usymbol(unit(q)) != DEFAULT_USYMBOL #Print the unit symbol if it exists (not default)
        if pretty
            return _pretty_print_quantity(io, q)
        else
            return _parsable_print_quantity(io, q)
        end 
    else
        return Base.show(io, ubase(q), pretty=pretty) #Print SI version (emphesizes that units are treated as SI quantities)
    end
end

#Fallback for showing pure units (don't pretty print)
function Base.show(io::IO, u::AbstractUnits; pretty=PRETTY_DIM_OUTPUT[])
    return Base.show_default(io, u)
end
    
function Base.show(io::IO, d::D; pretty=PRETTY_DIM_OUTPUT[]) where D<:AbstractDimensions{<:Real}
    dimnames = collect(static_fieldnames(D))
    dimunits = unit_symbols(D)

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

function Base.show(io::IO, d::D; pretty=PRETTY_DIM_OUTPUT[]) where D<: AbstractDimensions{<:Union{<:Real,Missing}}
    return show(io, missing)
end

function Base.show(io::IO, d::D; pretty=PRETTY_DIM_OUTPUT[]) where D<: AbstractDimensions
    vals = map(fn->fn=>getproperty(d,fn), dimension_names(D))
    return show(io, vals)
end

function _parsable_print_quantity(io::IO, q::AbstractQuantity)
    print(io, "(", ustrip(q), ")")
    return _parsable_print_unit(io, unit(q))
end

function _pretty_print_quantity(io::IO, q::AbstractQuantity)
    print(io, ustrip(q), " ")
    return _print_pretty_unit(io, unit(q))
end

function _parsable_print_unit(io::IO, u::AbstractUnitLike)
    usymb = usymbol(u)
    if usymb == DEFAULT_USYMBOL
        return Base.show(io, u, pretty=false)
    else
        return Base.print(io, usymb)
    end
end

function _print_pretty_unit(io::IO, u::AbstractUnitLike)
    usymb = usymbol(u)
    if usymb == DEFAULT_USYMBOL
        return Base.show(io, u, pretty=true)
    else
        return print(io, usymb)
    end
end


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

