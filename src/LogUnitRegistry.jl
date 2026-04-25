
module LogUnitRegistry

#RegistryTools contains all you need to build a registry in one simple import
using ..RegistryTools
import ..RegistryTools.dimtype 
import ..RegistryTools.unittype

const UNIT_LOCK = ReentrantLock()
const AffLogUnits = Union{Units{Dimensions{FixRat32}, AffineTransform{Float64}}, Units{Dimensions{FixRat32},ExpAffTransform{Float64}}}
const UNITS = PermanentDict{Symbol, AffLogUnits}()

#Additional conversion functions
Base.convert(::Type{AffLogUnits}, u::AbstractUnitLike) = convert(Units{Dimensions{FixRat32}, AffineTransform{Float64}}, u)
Base.convert(::Type{AffLogUnits}, u::Units{<:Any, <:ExpAffTransform}) = convert(Units{Dimensions{FixRat32}, ExpAffTransform{Float64}}, u)
Base.convert(::Type{AffLogUnits}, u::QuantUnion) = convert(Units{Dimensions{FixRat32}, AffineTransform{Float64}}, u)

#Fill the UNITS registry with default values
registry_defaults!(UNITS)

#Ueses a ReentrantLock() on register_unit to prevent race conditions when multithreading
register_unit(p::Pair{String, <:Any}) = lock(UNIT_LOCK) do 
    register_unit!(UNITS, p)
end

#Parsing functions that don't require a dictionary argument
uparse(str::String) = RegistryTools.uparse(str, UNITS)
qparse(str::String) = RegistryTools.qparse(str, UNITS)

#String macros are possible now that we are internally referring to UNITS
macro u_str(str)
    return suparse_expr(str, UNITS)
end

macro ud_str(str)
    return uparse_expr(str, UNITS)
end

macro q_str(str)
    return qparse_expr(str, UNITS)
end

macro U_str(str)
    suexpr = suparse_expr(str, UNITS)
    return :($typeof($suexpr))
end

macro D_str(str)
    suexpr = suparse_expr(str, UNITS)
    return :($dimtype($suexpr))
end

#Add these functions to facilitate knowing types ahead of time, DO NOT EXPORT IF MULTIPLE REGISTRIES ARE USED
utype() = RegistryTools.regunittype(UNITS)
dtype() = RegistryTools.regdimtype(UNITS)

#Registry is exported but these functions/macros are not (in case user wants their own verison)
#You can import these by invoking `using .Registry`
export @u_str, @ud_str, @q_str, @U_str, @D_str, uparse, qparse, register_unit

end