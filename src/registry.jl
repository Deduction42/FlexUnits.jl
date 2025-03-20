
module UnitRegistry

#RegistryTools contains all you need to build a registry in one simple import
using ..RegistryTools

const UNIT_LOCK = ReentrantLock()
const DIMS = DEFAULT_DIMENSONS
const UNITS = PermanentDict{Symbol, AffineUnits{DIMS}}()

#Fill the UNITS registry with default values
registry_defaults!(UNITS)

#Ueses a ReentrantLock() to prevent race conditions when multithreading
register_unit(p::Pair{String,<:AbstractUnitLike}) = lock(UNIT_LOCK) do 
    register_unit!(UNITS, p)
end

#Parsing functions that don't require adding a dictionary
uparse(str::String) = RegistryTools.uparse(str, UNITS)

#String macros are possible now that we are internally referring to UNITS
macro u_str(str)
    return esc(uparse_expr(str, UNITS))
end

#Registry is exported but these functions/macros are not (in case user wants their own verison)
#You can import these by invoking `using .Registry`
export @u_str, @q_str, uparse, qparse, register_unit

end