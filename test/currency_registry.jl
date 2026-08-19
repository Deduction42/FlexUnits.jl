module CurrencyUnits

using FlexUnits
using FlexUnits.RegistryTools

@kwdef struct MoneyDimensions{P} <: AbstractDimensions{P}
    m  ::P = zero(FixRat32)
    kg ::P = zero(FixRat32)
    s  ::P = zero(FixRat32)
    A  ::P = zero(FixRat32)
    K  ::P = zero(FixRat32)
    cd ::P = zero(FixRat32)
    mol::P = zero(FixRat32)
    EUR::P = zero(FixRat32)
end
MoneyDimensions(args::Real...) = MoneyDimensions{FixRat32}(args...)

const UNITS = PermanentDict{Symbol,Units{MoneyDimensions{FixRat32},AffineTransform{Float64}}}()

registry_defaults!(UNITS)
register_unit!(UNITS, "EUR" => MoneyDimensions{FixRat32}(EUR=1))
register_unit!(UNITS, "kEUR" => 1000*UNITS[:EUR])

uparse(str::String) = RegistryTools.uparse(str, UNITS)
qparse(str::String) = RegistryTools.qparse(str, UNITS)

macro u_str(str); return suparse_expr(str, UNITS); end
macro ud_str(str); return uparse_expr(str, UNITS); end
macro q_str(str); return qparse_expr(str, UNITS); end
macro U_str(str); suexpr = suparse_expr(str, UNITS); return :($typeof($suexpr)); end
macro D_str(str); suexpr = suparse_expr(str, UNITS); return :($dimtype($suexpr)); end

utype() = RegistryTools.regunittype(UNITS)
dtype() = RegistryTools.regdimtype(UNITS)

export @u_str, @ud_str, @q_str, @U_str, @D_str, uparse, qparse, utype, dtype

end