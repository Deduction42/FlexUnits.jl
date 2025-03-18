
module RegistryTools

import ..AbstractUnitLike, ..AbstractUnits, ..AbstractAffineUnits, ..AbstractScalarUnits, ..AbstractDimensions
import ..AffineUnits, ..ScalarUnits, ..Dimensions, ..RealQuantity, ..UnionQuantity
import ..uscale, ..uoffset, ..dimension, ..usymbol, ..asunit, ..constructorof

export 
    AbstractUnitLike, AbstractUnits, AbstractAffineUnits, AbstractScalarUnits, AbstractDimensions,
    AffineUnits, ScalarUnits, Dimensions, RealQuantity, UnionQuantity, uscale, uoffset, dimension, usymbol,
    PermanentDict, register_unit!, registry_defaults!, uparse, usparse, qparse, uparse_expr


"""
    PermanentDict{K,T}

A dict where any additions made are permanent and not allowed to be changed.
This prevents astonishing behaviour from string macros, which only look up the units
at compile time. If a value was changed since the lookup, the string macro will always
return the value looked up at compile time. Attempts to delete entries, or modify
them will result in an error. Reassinging the same value will result in no action.
"""
struct PermanentDict{K,V} <: AbstractDict{K,V}
    data :: Dict{K,V}
end
PermanentDict{K,V}() where{K,V} = PermanentDict{K,V}(Dict{K,V}())
Base.iterate(d::PermanentDict) = iterate(d.data)
Base.iterate(d::PermanentDict, n::Integer) = iterate(d.data, n)
Base.length(d::PermanentDict)  = length(d.data)

"""
    PermanentDictError

An error class that specializes in illegal operations on PermanentDict
"""
struct PermanentDictError <: Exception
    msg::String
end
Base.showerror(io::IO, e::PermanentDictError) = print(io, "PermanentDictError: ", e.msg)


Base.get(d::PermanentDict, k, v) = get(d.data, k, v)
Base.getindex(d::PermanentDict, k) = d.data[k]
function Base.setindex!(d::PermanentDict, v, k)
    if haskey(d.data, k)
        d.data[k] == v || throw(PermanentDictError("Key $(k) already exists with value $(d.data[k]), changing it to $(v) not allowed)"))
    else 
        d.data[k] = v 
    end
    return d
end
Base.delete!(d, k) = throw(PermanentDictError("Removing entries is prohibited"))

"""
    register_unit!(reg::AbstractDict{Symbol,<:AffineUnits}, p::Pair{String,<:AbstractUnitLike})

Registers a unit into the dictionary. This version should be used instead of the "Symbol" versions as 
it will validate whether or not the proposed name can be parsed.
"""
function register_unit!(reg::AbstractDict{Symbol,<:AffineUnits}, p::Pair{String,<:AbstractUnitLike})
    (k,v) = p
    name = Meta.parse(k)

    if name isa Symbol
        return _register_unit!(reg, name => v)
    else
        throw(ArgumentError("Invalid unit registration name for $(k), it did not parse as a Symbol"))
    end
end


function _register_unit!(reg::AbstractDict{Symbol,<:AffineUnits}, p::Pair{Symbol,<:AbstractDimensions})
    (k,v) = p
    return setindex!(reg, AffineUnits(dims=v, symbol=k), k)
end

function _register_unit!(reg::AbstractDict{Symbol,<:AffineUnits}, p::Pair{Symbol,<:AbstractUnitLike})
    (k,v) = p
    vn = AffineUnits(scale=uscale(v), offset=uoffset(v), dims=dimension(v), symbol=k)
    return setindex!(reg, vn, k)
end

function add_prefixes!(reg::AbstractDict{Symbol,<:AffineUnits{D}}, u::Symbol, prefixes::NamedTuple) where D<:AbstractDimensions
    original = reg[u]
    for (name, scale) in pairs(prefixes)
        newname = Symbol(string(name)*string(u))
        reg[newname] = AffineUnits{D}(scale=uscale(original)*scale, dims=dimension(original), symbol=newname)
    end
    return reg
end


function registry_defaults!(reg::AbstractDict{Symbol, AffineUnits{Dims}}) where Dims <: AbstractDimensions
    #reg = PermanentDict{Symbol, AffineUnits{DEFAULT_DIMENSONS}}()
    si_prefixes = (f=1e-15, p=1e-12, n=1e-9, μ=1e-6, u=1e-6, m=1e-3, c=1e-2, d=0.1, k=1e3, M=1e6, G=1e9, T=1e12)

    #SI dimensional units
    _register_unit!(reg, :m => Dims(length=1))
    _register_unit!(reg, :g => asunit(0.001*Dims(mass=1)))
    _register_unit!(reg, :t => asunit(1000*Dims(mass=1)))
    _register_unit!(reg, :s => Dims(time=1))
    _register_unit!(reg, :A => Dims(current=1))
    _register_unit!(reg, :K => Dims(temperature=1))
    _register_unit!(reg, :cd => Dims(luminosity=1))
    _register_unit!(reg, :mol => Dims(amount=1))
    
    add_prefixes!(reg, :m, si_prefixes[( :f, :p, :n, :μ, :u, :m, :c, :d, :k, :M, :G )])
    add_prefixes!(reg, :g, si_prefixes[( :n, :μ, :u, :m, :k)])
    add_prefixes!(reg, :t, si_prefixes[( :k, :M, :G)])
    add_prefixes!(reg, :s, si_prefixes[( :f, :p, :n, :μ, :m)])
    add_prefixes!(reg, :A, si_prefixes[( :n, :μ, :u, :m, :k)])
    add_prefixes!(reg, :K, si_prefixes[( :m, )])
    add_prefixes!(reg, :cd, si_prefixes[( :m, )])
    add_prefixes!(reg, :mol, si_prefixes[( :p, :n, :μ, :u, :m, :k )])

    #SI derived units
    m = reg[:m]
    kg = reg[:kg]
    s = reg[:s]
    A = reg[:A]
    mol = reg[:mol]
    K = reg[:K]

    _register_unit!(reg, :L => reg[:dm]^3)
    _register_unit!(reg, :Hz => inv(s))
    _register_unit!(reg, :N => kg*m/s^2); N = reg[:N]
    _register_unit!(reg, :Pa => N/m^2);
    _register_unit!(reg, :J => N*m); J = reg[:J]
    _register_unit!(reg, :W => J/s); W = reg[:W]
    _register_unit!(reg, :C => A*s); C = reg[:C]
    _register_unit!(reg, :V => W/A); V = reg[:V]
    _register_unit!(reg, :Ω => V/A); Ω = reg[:Ω]
    _register_unit!(reg, :F => C/V) 
    _register_unit!(reg, :ohm => Ω) 
    _register_unit!(reg, :S => A/V)
    _register_unit!(reg, :H => N*m/A^2)
    _register_unit!(reg, :T => N/(A*m))
    _register_unit!(reg, :Wb => V*s)

    _register_unit!(reg, :°C => AffineUnits(offset=273.15, dims=dimension(K)))    

end

function uparse(str::String, reg::AbstractDict{Symbol, <:AbstractUnitLike})
    return eval(uparse_expr(str, reg))
end

function usparse(str::String, reg::AbstractDict{Symbol,U}) where {U<:AbstractUnits}
    u = uparse(str, reg)
    return convert(ScalarUnits{dimtype(U)}, u)
end

function qparse(str::String, reg::AbstractDict{Symbol,U}) where {U<:AbstractUnits}
    u = uparse(str, reg)
    return convert(Quantity{Float64, dimtype(U)}, u)
end

function uparse_expr(str::String, reg::AbstractDict{Symbol, U}) where U <: AbstractUnitLike
    ex = uparse_expr(Meta.parse(str), reg)
    return :($change_symbol($(ex), $(str)))
end

function uparse_expr(ex::Expr, reg::AbstractDict{Symbol, U}) where U <: AbstractUnitLike
    if ex.head != :call
        throw(ArgumentError("Unexpected expression: $ex. Only `:call` is expected."))
    end
    ex.args[2:end] = map(Base.Fix1(lookup_unit_expr, reg), ex.args[2:end])

    #Make sure the return type is the same as the registry
    return :($convert($U, $ex))
end

function uparse_expr(ex::Symbol, reg::AbstractDict{Symbol, U}) where U <: AbstractUnitLike
    return lookup_unit_expr(reg, ex)
end

function lookup_unit_expr(reg::AbstractDict{Symbol,<:AbstractUnitLike}, ex::Expr)
    return uparse_expr(ex, reg)
end

function lookup_unit_expr(reg::AbstractDict{Symbol,<:AbstractUnitLike}, name::Symbol)
    if !haskey(reg, name)
        throw(ArgumentError("Symbol $sym not found in the provided unit registry."))
    else
        return reg[name]
    end
end

lookup_unit_expr(reg::AbstractDict{Symbol,<:AbstractUnitLike}, ex::Any) = ex

change_symbol(u::AbstractUnitLike, s::String) = change_symbol(u, Symbol(s))

function change_symbol(u::U, s::Symbol) where U<:AbstractAffineUnits
    return constructorof(U)(scale=uscale(u), offset=uoffset(u), dims=dimension(u), symbol=s)
end

function change_symbol(u::U, s::Symbol) where U<:AbstractScalarUnits
    return constructorof(U)(scale=uscale(u), dims=dimension(u), symbol=s)
end

end


#=================================================================================================
Design Decisions:
- Registry-agnostic parsing functions live in RegistryTools module 
    uparse(str, reg::Dict) parses `str` and fills in units according to `reg`, returning units
        variants will convert to appropriate type
    uparse_expr(str, reg::Dict) parses `str` like {x} and returns an expression
        variants add the conversion expression
    macros rely on the module environment for `reg`
        makes uparse_expr calls underneath

- RegistryTools should import all neccessary functions/types from the environment, and exports them
    using ..RegistryTools
    This simplifies the task of creating a new registry

- Registry-internal parsing functions live inside the Registry module (they are not automatically exported)
    uparse(str::String) (and variants)
    @u_str (and variaents)
    register_unit(p::Pair{Symbol, <:AbstractUnitLike})

- Parsing will have three variants (not automatically exported, )
    uparse will return whatever units are in the dictionary (AffineUnits{D})
        u"..." is the equivalent macro
    suparse will convert uparse results into ScalarUnits{D}
        su"..." is the 
    qparse will convert uparse results into RealQuantity{Float64, D}
        q"..." is the equivalent macro

- Registry-internal registering function will need to have a ReentrantLock applied

=================================================================================================#

#=
#Test code
reg = RegistryTools.PermanentDict{Symbol, AffineUnits{DEFAULT_DIMENSONS}}()
RegistryTools.registry_defaults!(reg)
5*RegistryTools.uparse("°C", reg)
=#
