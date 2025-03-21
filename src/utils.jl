

function Base.show(io::IO, q::UnionQuantity{<:Any, <:AbstractDimensions}; pretty=PRETTY_DIM_OUTPUT[])
    if pretty
        return _print_pretty_space(io, q)
    else
        return _print_string_macro(io, q)
    end
end

function Base.show(io::IO, q::UnionQuantity{<:Any, <:AbstractUnits}; pretty=PRETTY_DIM_OUTPUT[])
    if usymbol(unit(q)) != DEFAULT_USYMBOL #Print the unit symbol if it exists (not default)
        if pretty
            return _print_pretty_space(io, q)
        else
            return _print_string_macro(io, q)
        end 
    else
        return Base.show(io, ubase(q), pretty=pretty) #Print SI version (emphesizes that units are treated as SI quantities)
    end
end

#Fallback for showing pure units (don't pretty print)
function Base.show(io::IO, u::AbstractUnits; pretty=PRETTY_DIM_OUTPUT[])
    return Base.show_default(io, u)
end
    
function Base.show(io::IO, d::D; pretty=PRETTY_DIM_OUTPUT[]) where D<: AbstractDimensions
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

    num = filter(k-> d[k]>0, dimnames)
    den = filter(k-> d[k]<0, dimnames)

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

function _print_string_macro(io::IO, q::UnionQuantity)
    print(io, "(", ustrip(q), ")")
    return _print_unit_macro(io, unit(q))
end

function _print_pretty_space(io::IO, q::UnionQuantity)
    print(io, ustrip(q), " ")
    return _print_pretty_unit(io, unit(q))
end

function _print_unit_macro(io::IO, u::AbstractUnitLike)
    usymb = usymbol(u)
    print(io, "u\"")

    if usymb == DEFAULT_USYMBOL
        Base.show(io, u, pretty=false)
    else
        Base.print(io, usymb)
    end
    return print(io, "\"")
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

