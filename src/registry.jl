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

const PREFIXES = (f=1e-15, p=1e-12, n=1e-9, μ=1e-6, u=1e-3, c=1e-2, d=0.1, k=1e3, M=1e6, G=1e9, T=1e12)

function register_unit!(reg::AbstractDict{Symbol,<:AffineUnits}, p::Pair{Symbol,<:AbstractDimensions})
    (k,v) = p
    return setindex!(reg, AffineUnits(dims=v, symbol=k), k)
end

function register_unit!(reg::AbstractDict{Symbol,<:AffineUnits}, p::Pair{Symbol,<:AbstractUnitLike})
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



function fill_default_units!(reg::AbstractDict{AffineUnits{D}}) where D <: AbstractDimensions
    reg = PermanentDict{Symbol, AffineUnits{DEFAULT_DIMENSONS}}()
    register_unit!(reg, :m => Dimensions(mass=1))
    register_unit!(reg, :g => ScalarUnits(scale=0.001, dims=Dimensions(mass=1)))
    register_unit!(reg, :s => Dimensions(time=1))
    register_unit!(reg, :A => Dimensions(current=1))
    register_unit!(reg, :K => Dimensions(temperature=1))
    register_unit!(reg, :cd => Dimensions(luminosity=1))
    register_unit!(reg, :mol => Dimensions(amount=1))

    add_prefixes!(reg, :m, PREFIXES[( :f, :p, :n, :μ, :u, :c, :d, :k, :M, :G )])
end


