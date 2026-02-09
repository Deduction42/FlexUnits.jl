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

function pretty_print_units(x::Bool) 
    PRETTY_DIM_OUTPUT[] = x 
    return x 
end

#Dispatching show methods
Base.show(io::IO, q::QuantUnion) = qshow(io, q, pretty=PRETTY_DIM_OUTPUT[])
Base.string(q::QuantUnion) = qstring(q, pretty=PRETTY_DIM_OUTPUT[])
Base.show(io::IO, u::AbstractUnitLike) = ushow(io, u, pretty=PRETTY_DIM_OUTPUT[])
Base.string(u::AbstractUnitLike) = ustring(u, pretty=PRETTY_DIM_OUTPUT[])

#=============================================================================================
Displaying quantities
=============================================================================================#
function qstring(x::QuantUnion; pretty=PRETTY_DIM_OUTPUT[])
    tmp_io = IOBuffer()
    qshow(tmp_io, x, pretty=pretty)
    return String(take!(tmp_io))
end

#Generic fallback for showing any uncaught quantity type
qshow(io::IO, q::QuantUnion; pretty=PRETTY_DIM_OUTPUT[]) = Base.show_default(io, q)

function qshow(io::IO, q::QuantUnion{<:Any, <:AbstractUnits}; pretty=PRETTY_DIM_OUTPUT[])
    if usymbol(unit(q)) != DEFAULT_USYMBOL #Print the unit symbol if it exists (not default)
        if pretty
            return qshow_pretty(io, q)
        else
            return qshow_parsable(io, q)
        end 
    else
        return qshow(io, ubase(q), pretty=pretty) #Print SI version (emphesizes that units are treated as SI quantities)
    end
end

function qshow(io::IO, q::QuantUnion{<:Any, <:AbstractDimLike}; pretty=PRETTY_DIM_OUTPUT[])
    if pretty
        return qshow_pretty(io, q)
    else
        return qshow_parsable(io, q)
    end
end

function qshow_parsable(io::IO, q::QuantUnion)
    print(io, "(", ustrip(q), ")")
    return ushow(io, unit(q), pretty=false)
end

function qshow_pretty(io::IO, q::QuantUnion)
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


#=
function ushow(io::IO, d::D; pretty=PRETTY_DIM_OUTPUT[]) where D<: AbstractDimensions{<:Union{<:Real,Missing}}
    return show(io, missing)
end
=#

#=
function ushow(io::IO, u::Units; pretty=PRETTY_DIM_OUTPUT[])
    if usymbol(u) != DEFAULT_USYMBOL
        return print(io, usymbol(u))
    else
        print(io, "Units(todims=$(todims(u)), dims=")
        show(io, dimension(u); pretty)
        return print(io, ")")
    end
end
=#

#=
function ushow(io::IO, ::Type{MirrorDims{D}}; pretty=PRETTY_DIM_OUTPUT[]) where {D<:AbstractDimensions}
    return print(io, "MirrorDims{$(D)}")
end
=#
#=
function ushow(io::IO, ::Type{Union{D,MirrorDims{D}}}; pretty=PRETTY_DIM_OUTPUT[]) where {D<:AbstractDimensions}
    return print(io, "MirrorUnion{$(D)}")
end
=#