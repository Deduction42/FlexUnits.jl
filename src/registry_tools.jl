
module RegistryTools

import ..AbstractUnitLike, ..AbstractUnits, ..AbstractAffineUnits, ..AbstractDimensions
import ..AffineUnits, ..Dimensions, ..RealQuantity, ..UnionQuantity, ..DEFAULT_DIMENSONS
import ..uscale, ..uoffset, ..dimension, ..usymbol, ..asunit,  ..ubase, ..constructorof, ..dimtype

export 
    AbstractUnitLike, AbstractUnits, AbstractAffineUnits, AbstractDimensions, DEFAULT_DIMENSONS,
    AffineUnits, Dimensions, RealQuantity, UnionQuantity, uscale, uoffset, dimension, usymbol,
    PermanentDict, register_unit!, registry_defaults!, uparse, qparse, uparse_expr, qparse_expr, dimtype


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

function _register_unit!(reg::AbstractDict{Symbol,<:AffineUnits}, p::Pair{Symbol,<:AbstractUnitLike})
    (k,v) = p
    vn = AffineUnits(scale=uscale(v), offset=uoffset(v), dims=dimension(v), symbol=k)
    return setindex!(reg, vn, k)
end

function _register_unit!(reg::AbstractDict{Symbol,<:AffineUnits}, p::Pair{Symbol, <:UnionQuantity})
    return _register_unit!(reg, p[1]=>asunit(p[2]))
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
    
    _register_unit(p::Pair) = _register_unit!(reg, p)
    add_prefixes(symb::Symbol, prfx::NamedTuple) = add_prefixes!(reg, symb, prfx)

    #SI dimensional units
    _register_unit(:NoDims => Dims())
    _register_unit(:m => Dims(length=1))
    _register_unit(:g => 0.001*Dims(mass=1))
    _register_unit(:t => 1000*Dims(mass=1))
    _register_unit(:s => Dims(time=1))
    _register_unit(:A => Dims(current=1))
    _register_unit(:K => Dims(temperature=1))
    _register_unit(:cd => Dims(luminosity=1))
    _register_unit(:mol => Dims(amount=1))
    
    add_prefixes(:m, si_prefixes[( :f, :p, :n, :μ, :u, :m, :c, :d, :k, :M, :G )])
    add_prefixes(:g, si_prefixes[( :n, :μ, :u, :m, :k)])
    add_prefixes(:t, si_prefixes[( :k, :M, :G)])
    add_prefixes(:s, si_prefixes[( :f, :p, :n, :μ, :m)])
    add_prefixes(:A, si_prefixes[( :n, :μ, :u, :m, :k)])
    add_prefixes(:K, si_prefixes[( :m, )])
    add_prefixes(:cd, si_prefixes[( :m, )])
    add_prefixes(:mol, si_prefixes[( :μ, :u, :m, :k )])

    #SI derived units
    m = dimension(reg[:m])
    kg = dimension(reg[:kg])
    s = dimension(reg[:s])
    A = dimension(reg[:A])
    mol = dimension(reg[:mol])
    K = dimension(reg[:K])

    _register_unit(:percent => 0.01*reg[:NoDims])
    _register_unit(:L => reg[:dm]^3)
    _register_unit(:Hz => inv(s))
    _register_unit(:N => kg*m/s^2); N = reg[:N]
    _register_unit(:Pa => N/m^2);
    _register_unit(:J => N*m); J = reg[:J]
    _register_unit(:W => J/s); W = reg[:W]
    _register_unit(:C => A*s); C = reg[:C]
    _register_unit(:V => W/A); V = reg[:V]
    _register_unit(:Ω => V/A); Ω = reg[:Ω]
    _register_unit(:F => C/V) 
    _register_unit(:ohm => Ω) 
    _register_unit(:S => A/V)
    _register_unit(:H => N*m/A^2)
    _register_unit(:T => N/(A*m))
    _register_unit(:Wb => V*s)

    add_prefixes(:L, si_prefixes[(:μ, :u, :m)])
    add_prefixes(:Hz, si_prefixes[(:n, :μ, :u, :m, :k, :M, :G)])
    add_prefixes(:N, si_prefixes[(:μ, :u, :m, :k)])
    add_prefixes(:Pa, si_prefixes[(:k,:M, :G)])
    add_prefixes(:J, si_prefixes[(:k,:M,:G)])
    add_prefixes(:W, si_prefixes[(:m, :k, :M, :G)])
    add_prefixes(:V, si_prefixes[(:p, :n, :μ, :u, :m, :k, :M, :G)])
    add_prefixes(:F, si_prefixes[(:f, :p, :n, :μ, :u, :m)])
    add_prefixes(:Ω, si_prefixes[(:n, :μ, :u, :m, :k, :M, :G)])
    add_prefixes(:ohm, si_prefixes[(:n, :μ, :u, :m, :k, :M, :G)])
    add_prefixes(:S, si_prefixes[(:n, :μ, :u, :m, :k, :M, :G)])
    add_prefixes(:Wb, si_prefixes[(:n, :μ, :u, :m)])

    #Common time units
    _register_unit(:min => 60*s); minute=reg[:min]
    _register_unit(:minute => minute)
    _register_unit(:h => 60*minute); h = reg[:h]
    _register_unit(:hr => h);
    _register_unit(:day => 24*h); day = reg[:day]
    _register_unit(:d => day); 
    _register_unit(:wk => 7*day);
    _register_unit(:yr => 365.25*day);

    #Common imperial units
    _register_unit(:inch => 2.54*reg[:cm]); inch = reg[:inch]
    _register_unit(:ft => 12*inch); ft = reg[:ft]
    _register_unit(:mi => 5280*ft)
    _register_unit(:lb => 0.453592*kg); lb = reg[:lb]
    _register_unit(:oz => (1/16)*lb)
    _register_unit(:psi => lb/inch^2)
    _register_unit(:lbf => 4.44822*N)
    _register_unit(:fl_oz => 29.5735*reg[:mL])
    _register_unit(:cup => 8*reg[:fl_oz])
    _register_unit(:pint => 2*reg[:cup])
    _register_unit(:quart => 2*reg[:pint])

    #Strictly affine temperature measurements
    _register_unit(:°C => AffineUnits(offset=273.15, dims=K))    
    _register_unit(:°F => AffineUnits(scale=5/9, offset=(273.15 - 32*5/9), dims=K))
    _register_unit(:degC => reg[:°C])
    _register_unit(:degF => reg[:°F])

    #Angled units: For the SI system, considers "degrees/radians" to be dimensionless 
    _register_unit(:rad => reg[:NoDims])
    _register_unit(:deg => (2*π/360)*reg[:NoDims])
    _register_unit(:rpm => (2*π/60)*(reg[:Hz]))

    #If you want to add "angle" as a dimension, you can overload the appropriate function

    #function apply_trig_func(f, q::UnionQuantity{<:Any, <:RadDimensions{T}}) where T
    #   baseq = ubase(q)
    #   assert_radians(unit(baseq))
    #   return quantity(f(ustrip(baseq)), RadDimensions{T}())
    #end

    #sin(q::UnionQuantity{<:RadDimensions}) = apply_trig_func(sin, q)

    return reg
end

function uparse(str::String, reg::AbstractDict{Symbol, <:AbstractUnitLike})
    return eval(uparse_expr(str, reg))
end

function qparse(str::String, reg::AbstractDict{Symbol,U}) where {U<:AbstractUnits}
    return eval(qparse_expr(str, reg))
end

function uparse_expr(str::String, reg::AbstractDict{Symbol, U}) where U <: AbstractUnitLike
    parsed = Meta.parse(str)
    if !(parsed isa Union{Expr,Symbol})
        throw(ArgumentError("Unexpected expression: String input \"$(str)\" was not parsed as an Expr or Symbol"))
    else
        ex = uparse_expr(parsed, reg)
        return :($change_symbol($(ex), $(str)))
    end
end

function qparse_expr(str::String, reg::AbstractDict{Symbol, U}) where U <: AbstractUnitLike
    parsed = Meta.parse(str)
    if !(parsed isa Union{Expr,Symbol})
        throw(ArgumentError("Unexpected expression: String input \"$(str)\" was not parsed as an Expr or Symbol"))
    else
        return qparse_expr(parsed, reg)
    end
end

function uparse_expr(ex::Union{Expr,Symbol}, reg::AbstractDict{Symbol, U}) where U <: AbstractUnitLike
    ex = _parse_expr(ex, reg)
    return :($convert($U, $ex))
end

function qparse_expr(ex::Union{Expr,Symbol}, reg::AbstractDict{Symbol, U}) where U
    ex = _parse_expr(ex, reg)
    return :(ubase($ex))
end

function _parse_expr(ex::Expr, reg::AbstractDict{Symbol, U}) where U <: AbstractUnitLike
    if ex.head != :call
        throw(ArgumentError("Unexpected expression: $ex. Only `:call` is expected."))
    end
    ex.args[2:end] = map(Base.Fix1(lookup_unit_expr, reg), ex.args[2:end])
    return ex
end

function _parse_expr(ex::Symbol, reg::AbstractDict{Symbol, U}) where U <: AbstractUnitLike
    return lookup_unit_expr(reg, ex)
end

function lookup_unit_expr(reg::AbstractDict{Symbol,<:AbstractUnitLike}, ex::Expr)
    return _parse_expr(ex, reg)
end

function lookup_unit_expr(reg::AbstractDict{Symbol,<:AbstractUnitLike}, name::Symbol)
    if !haskey(reg, name)
        throw(ArgumentError("Symbol $(name) not found in the provided unit registry."))
    else
        return reg[name]
    end
end

lookup_unit_expr(reg::AbstractDict{Symbol,<:AbstractUnitLike}, ex::Any) = ex

change_symbol(u::AbstractUnitLike, s::String) = change_symbol(u, Symbol(s))

function change_symbol(u::U, s::Symbol) where U<:AbstractAffineUnits
    return constructorof(U)(scale=uscale(u), offset=uoffset(u), dims=dimension(u), symbol=s)
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
