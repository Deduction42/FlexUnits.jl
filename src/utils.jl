function Base.show(io::IO, d::D) where D<: AbstractDimensions
    dimnames = collect(static_fieldnames(D))
    dimunits = unit_symbols(D)

    abs_dim_pwr(dim::Symbol) = pretty_unit_pwr(dimunits[dim], abs(d[dim]))

    function print_dim_series(series_io::IO, dims)
        if length(dims) == 1
            print(series_io, abs_dim_pwr(dims[1]))
        elseif length(dims) > 1
            series_str = "(" * join(abs_dim_pwr.(dims), "*") * ")"
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

function Base.show(io::IO, u::U) where U<:AbstractScalarUnits
    if usymbol(u) != DEFAULT_USYMBOL
        print(io, usymbol(u))
    else
        ustring = "$(uscale(u))*{$(dimension(u)))"
        print(io, ustring)
    end

    return nothing
end

function Base.show(io::IO, u::U) where U<:AbstractAffineUnits
    if usymbol(u) != DEFAULT_USYMBOL
        print(io, usymbol(u))
    else
        offset = uoffset(u)
        sgn = ifelse(offset >= 0, '+', '-')
        ustring = "$(uscale(u))*{$(dimension(u)) $(sgn) $(offset)}"
        print(io, ustring)
    end

    return nothing
end

function Base.show(io::IO, q::UnionQuantity{<:Any, <:AbstractDimensions})
    return print(io, ustrip(q), " ", unit(q))
end

function Base.show(io::IO, q::UnionQuantity{<:Any, <:AbstractUnits})
    symb = usymbol(unit(q))
    if symb != DEFAULT_USYMBOL
        print(io, ustrip(q), " ", unit(q))
    else
        print(io, ubase(q)) #Print SI version (emphesizes that units are treated as SI quantities)
    end

    return nothing
end



function pretty_unit_pwr(u::Symbol, p::Real)
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