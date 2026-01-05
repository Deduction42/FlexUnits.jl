module FlexUnitfulExt

import FlexUnits
import FlexUnits.UnitRegistry
import Unitful
import Unitful: @u_str


@generated function validate_upreferred()
    si_units = (length=u"m", mass=u"kg", time=u"s", current=u"A", temperature=u"K", luminosity=u"cd", amount=u"mol")

    for k in keys(si_units)
        Unitful.upreferred(si_units[k]) == si_units[k] || error("Found custom `Unitful.preferunits`: FlexUnits only supporsts the default SI option for `Unitful.upreferred`")
    end
    
    return true
end

function FlexUnits.Quantity(q::Unitful.Quantity)
    validate_upreferred()
    flex_dim   = convert(UnitRegistry.dimtype(), Unitful.dimension(q))
    base_scale = Unitful.ustrip(Unitful.upreferred(q))
    return FlexUnits.Quantity(base_scale, flex_dim)
end

function FlexUnits.Quantity{T}(q::Unitful.Quantity) where {T}
    flexq = FlexUnits.Quantity(q)
    return FlexUnits.Quantity{T}(FlexUnits.ustrip(flexq), FlexUnits.unit(flexq))
end

function FlexUnits.Quantity{T,U}(q::Unitful.Quantity) where {T,U}
    flexq = FlexUnits.Quantity(q)
    return FlexUnits.Quantity{T,U}(FlexUnits.ustrip(flexq), FlexUnits.unit(flexq))
end

Base.convert(::Type{FlexUnits.Quantity{T,U}}, q::Unitful.Quantity) where {T,U<:FlexUnits.Units} = FlexUnits.Quantity{T,U}(q)
Base.convert(::Type{FlexUnits.Quantity{T,U}}, q::Unitful.Quantity) where {T,U<:FlexUnits.AbstractDimensions} = FlexUnits.Quantity{T,U}(q)

function FlexUnits.uconvert(u::FlexUnits.AbstractUnitLike, q::Unitful.Quantity)
    check_dimensions(Unitful.unit(q), FlexUnits.dimension(u))
    base_unitful = Unitful.upreferred(q)
    base_flex = Unitful.ustrip(base_unitful)*FlexUnits.dimension(u)
    return FlexUnits.uconvert(u, base_flex)
end

function FlexUnits.uconvert(u::Unitful.Units, q::FlexUnits.Quantity)
    check_dimensions(u, FlexUnits.dimension(q))
    base_flex = FlexUnits.ubase(q)
    base_unitful = FlexUnits.ustrip(base_flex)*Unitful.upreferred(u)
    return Unitful.uconvert(u, base_unitful)
end

function check_dimensions(u::Unitful.Units, flexdims::FlexUnits.Dimensions)
    return check_dimensions(Unitful.dimension(u), flexdims)
end

function check_dimensions(dims::Unitful.Dimensions{D}, flexdims::FlexUnits.Dimensions{R}) where {D, R}
    validate_upreferred()
    if convert(FlexUnits.Dimensions{R}, dims) == flexdims
        return true
    else
        throw(FlexUnits.DimensionError((Unitful.Dimensions{D}(), flexdims)))
    end
end

function Base.convert(::Type{FlexUnits.Dimensions{R}}, d::Unitful.Dimensions{D}) where {R, D}
    FlexDim = FlexUnits.Dimensions{R}
    validate_upreferred()
    return prod(Base.Fix1(convert, FlexDim), D; init=FlexDim())
end

function Base.convert(::Type{FlexUnits.Dimensions{R}}, dims::Unitful.Dimension{D}) where {R, D}
    FlexDim = FlexUnits.Dimensions{R}
    validate_upreferred()

    D == :Length && return FlexDim(length=dims.power)
    D == :Mass && return FlexDim(mass=dims.power)
    D == :Time && return FlexDim(time=dims.power)
    D == :Current && return FlexDim(current=dims.power)
    D == :Temperature && return FlexDim(temperature=dims.power)
    D == :Luminosity && return FlexDim(luminosity=dims.power)
    D == :Amount && return FlexDim(amount=dims.power)
    error("Unknown dimension: $D")
end

end