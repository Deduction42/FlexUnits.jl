const DEFAULT_RATIONAL = FixRat32
const DEFAULT_USYMBOL = :_

#Arithmetic Types in Base that support +-*/ (use instead of Any to avoid ambiguity)
const MathUnion = Union{Real, Complex, AbstractArray}
const NumUnion = Union{Real, Complex}

"""
    AbstractUnitLike

A broad class representing anything that can be interpreted as a unit or dimension
"""
abstract type AbstractUnitLike end

"""
    AbstractUnitTransform

An abstract object representing a unit conversion formula. 
Any object that subtypes this is made callable.

```julia
# Callable form 
utrans = uconvert(u"°C", u"°F")
utrans(0.0)
31.999999999999986

# Shorthand callable form (syntactic sugar)
(u"°C" |> u"°F")(0.0)
31.999999999999986
```
"""
abstract type AbstractUnitTransform end

"""
    AbstractUnitLike

A broad class representing anything that can be interpreted as a unit
"""
abstract type AbstractDimLike <: AbstractUnitLike end

"""
    AbstractDimensions

A class that represents a specific dimensional schema, by default, its
only child is Dimensions (SI units), but users can build their own versions 
(even based on imperial measurements, or where angles are a dimension)
"""
abstract type AbstractDimensions{P} <: AbstractDimLike end

"""
    AbstractUnits

All units in FlexDims contain two entities: 
(1) AbstractDimLike (anything that can be interpreted as a dimension)
(2) AbstractTransform (contains the formula to convert unit value to dimensions)
"""
abstract type AbstractUnits{D, T<:AbstractUnitTransform} <: AbstractUnitLike end


const UnitOrDims{D} = Union{D, AbstractUnits{D}} where D<:AbstractDimensions
const ScalarOrVec{T} = Union{T, AbstractVector{T}} where T
abstract type AbstractUnitMap{U<:UnitOrDims{<:AbstractDimensions}} end

#Base.@assume_effects :consistent static_fieldnames(t::Type) = Base.fieldnames(t)
static_fieldnames(t::Type) = Base.fieldnames(t)
Base.eltype(::Type{<:AbstractDimensions{P}}) where P = P

#=======================================================================================
AbstractDimensions API
=======================================================================================#
#Dimension constructor from other types of dimensions
function (::Type{D})(x::AbstractDimensions) where {P, D<:AbstractDimensions{P}}
    return D(map(Base.Fix1(getproperty, x), static_fieldnames(D))...)
end

#Assign a single value to all dimensions
function (::Type{D})(x::Real) where {P, D<:AbstractDimensions{P}}
    f(anything) = x
    return D(map(f, static_fieldnames(D))...)
end

dimension(u::AbstractDimLike) = u
dimval(u::AbstractDimensions) = u
usymbol(u::AbstractDimLike) = DEFAULT_USYMBOL
uscale(u::AbstractUnitLike) = uscale(todims(u))
uoffset(u::AbstractUnitLike) = uoffset(todims(u))
dimtype(x::Any) = dimtype(typeof(x))
dimtype(::Type{T}) where T = error("dimtype not yet implemented for type $(T)")
dimtype(::Type{D}) where D<:AbstractDimLike = D
dimtype(::Type{U}) where {D, U<:AbstractUnits{D}} = dimtype(D)
dimvaltype(x::Any) = dimvaltype(typeof(x)) 
dimvaltype(::Type{T}) where T = dimtype(T)
dimpowtype(::Type{D}) where {P, D<:AbstractDimensions{P}} = P
Base.getindex(d::AbstractDimensions, k::Symbol) = getproperty(d, k)
dimpowtype(::Type{U}) where {U<:AbstractUnitLike} = dimpowtype(dimtype(U))
is_dimension(u::AbstractUnitLike) = is_identity(todims(u))
is_scalar(u::AbstractUnitLike) = is_scalar(todims(u))
udynamic(u::AbstractDimensions) = u

#=======================================================================================
Basic SI dimensions and transforms
=======================================================================================#
"""
NoTransform object, the default transform returned by todims(x::AbstractDimensionLike). Calling it results in 
an identity.
```julia
t = NoTransform()
t("anything")
"anything"
```
"""
struct NoTransform <: AbstractUnitTransform end 
Base.broadcastable(utrans::AbstractUnitTransform) = Ref(utrans)
(t2::AbstractUnitTransform)(t1::AbstractUnitTransform) = t2 ∘ t1
(t::NoTransform)(x) = x
(t::NoTransform)(t0::AbstractUnitTransform) = t0

Base.:∘(t1::NoTransform, t2::AbstractUnitTransform) = t2 
Base.:∘(t1::AbstractUnitTransform, t2::NoTransform) = t1 
Base.:∘(t1::NoTransform, t2::NoTransform) = t1
Base.inv(t::NoTransform) = t
uscale(t::NoTransform)  = 1
uoffset(t::NoTransform) = 0
todims(u::AbstractDimLike) = NoTransform()
is_identity(t::NoTransform) = true
is_scalar(t::NoTransform) = true

"""
    Dimensions{P}

Basic SI dimensions:
    length = m, 
    mass = kg, 
    time = s, 
    current = A, 
    temperature = K, 
    luminosity = cd, 
    amount = mol
"""
@kwdef struct Dimensions{P} <: AbstractDimensions{P}
    length::P = FixedRational(0)
    mass::P = FixedRational(0)
    time::P = FixedRational(0)
    current::P = FixedRational(0)
    temperature::P = FixedRational(0)
    luminosity::P = FixedRational(0)
    amount::P = FixedRational(0)
end
Dimensions(args...) = Dimensions{FixRat32}(args...)

function unit_symbols(::Type{<:Dimensions})
    return Dimensions{Symbol}(
        length=:m, mass=:kg, time=:s, current=:A, temperature=:K, luminosity=:cd, amount=:mol
    )
end

"""
    dimension_names(::Type{<:AbstractDimensions})

Overloadable fieldnames for a dimension (returns a tuple of the dimension's field names)
This should be static so that it can be hardcoded during compilation. 
Can use this to overload the default "fieldnames" behaviour
"""
@inline function dimension_names(::Type{D}) where {D<:AbstractDimensions}
    return static_fieldnames(D)
end

"""
    StaticDims{D} <: AbstractDimLike

Static dimensions where the `D`` is the dimension value. This improves performance when dimensions are
statically inferrable.
"""
struct StaticDims{D} <: AbstractDimLike
    function StaticDims{D}() where D
        return (D isa AbstractDimensions) ? new{D}() : error("Type parameter must be a dimension")
    end
end 
StaticDims(D::AbstractDimensions) = StaticDims{D}()
StaticDims{D}(d::AbstractDimensions) where D = (D == d) ? StaticDims{D} : throw(ArgumentError("Dimesion $(d) must be equal to $(D)"))
dimval(::Type{<:StaticDims{D}}) where D = D
dimval(d::StaticDims) = dimval(typeof(d))
udynamic(u::StaticDims{D}) where D = D
ustatic(u::AbstractDimensions) = StaticDims{u}()
dimvaltype(::Type{D}) where {d, D<:StaticDims{d}} = dimtype(d)

"""
    NoDims

Placehoder for a unitless dimension in cases where the base dimension type is uninferrable
"""
@kwdef struct NoDims <: AbstractDimensions{Bool}
    any :: Bool = false #Has a single value set to false so that it doesn't trigger "isunknown(NoDims())=true")
    NoDims(x::Bool) = new(x)
end
Base.getproperty(d::NoDims, ::Symbol) = getfield(d,:any)
dimension(x::NumUnion) = NoDims()
dimtype(::Type{<:NumUnion}) = NoDims
unit_symbols(::Type{<:NoDims}) = NoDims()
ustrip(x::NumUnion) = x 
dstrip(x::NumUnion) = x

#=======================================================================================
Affine Units and Transforms
=======================================================================================#
"""
    @kwdef struct AffineTransform{T<:Real} <: AbstractUnitTransform
        scale  :: T = 1.0
        offset :: T = 0.0
    end

A type representing an affine transfomration formula that can be
used to convert values from one affine unit to another. This object is callable.

# Constructors
    AffineTransform(scale::Real, offset::Real)
    AffineTransform(; scale, offset)
"""
@kwdef struct AffineTransform{T<:Real} <: AbstractUnitTransform
    scale  :: T = 1.0
    offset :: T = 0.0
end
AffineTransform(scale::T1, offset::T2) where {T1,T2} = AffineTransform{promote_type(T1,T2)}(scale, offset)
(t::AffineTransform)(x) = muladd(x, t.scale, t.offset)
(t::AffineTransform)(x::AbstractArray) = t.(x)
(t::AffineTransform)(x::Tuple) = map(t, x)
(t::AffineTransform)(t0::AbstractUnitTransform) = t ∘ t0

function Base.:∘(t2::AffineTransform, t1::AffineTransform)
    return AffineTransform(
        scale  = t1.scale*t2.scale,
        offset = t2.offset + t2.scale*t1.offset 
    )
end
Base.inv(t::AffineTransform) = AffineTransform(scale=inv(t.scale), offset=-t.offset/t.scale)

uscale(t::AffineTransform) = t.scale 
uoffset(t::AffineTransform) = t.offset
is_identity(t::AffineTransform) = isone(t.scale) & iszero(t.offset)
is_scalar(t::AffineTransform) = iszero(t.offset)
remove_offset(t::AffineTransform) = AffineTransform(scale=t.scale, offset=0)

"""
    @kwdef struct Units{D<:AbstractDimLike, T<:AbstractUnitTransform} <: AbstractUnits{D, T}
        dims   :: D
        todims :: T
        symbol :: Symbol = DEFAULT_USYMBOL
    end

A dynamic unit object that contains dimensions (dims) and its conversion formula to said dimensions (todims). The conversion 
formula determines what kind of unit is referred to. An AffineTransform implies affine units, a NoTransform implies dimensions.
Units with dynamic dimensions can generated through the `@ud_str` macro, while units with static dimensions can be generated
throiugh the `@u_str` macro.

# Constructors
    Units{D}(units, todims::T, symbol=DEFAULT_USYMBOL) where {D,T<:AbstractUnitTransform} 
    Units(dims::D, todims::AbstractUnitTransform=NoTransform(), symbol=DEFAULT_USYMBOL) where D<:AbstractDimensions 
    Units(units::D, todims::AbstractUnitTransform, symbol=DEFAULT_USYMBOL) where D<:AbstractUnits

```julia
julia> 1*(5ud"°C") #Operations on units eagerly convert to dimensions
278.15 K

julia> 1*(5du"°C") |> ud"°C" #Converts operation results back to Celsius
5.0 °C

julia> (5ud"°C" + 2ud"°C") |> ud"°C" #Operation adds values in Kelvin, results converted back to Celsius
280.15 °C

julia> (ustrip(5ud"°C") + ustrip(2ud"°C"))*u"°C" #Strips, adds raw quantity values, converts raw number to Celsius
7 °C
```
"""
@kwdef struct Units{D<:AbstractDimLike, T<:AbstractUnitTransform} <: AbstractUnits{D, T}
    dims   :: D
    todims :: T
    symbol :: Symbol = DEFAULT_USYMBOL
end
Units{D}(units, todims::T, symbol=DEFAULT_USYMBOL) where {D,T<:AbstractUnitTransform} = Units{D,T}(units, todims, symbol)
Units(dims::D, todims::AbstractUnitTransform=NoTransform(), symbol=DEFAULT_USYMBOL) where D<:AbstractDimLike = Units(dims, todims, symbol)
Units(units::U, todims::AbstractUnitTransform, symbol=DEFAULT_USYMBOL) where U<:AbstractUnits = Units(dimension(assert_dimension(units)), todims, symbol)

todims(u::Units) = u.todims
dimension(u::Units) = u.dims 
usymbol(u::Units) = u.symbol
remove_offset(u::U) where U<:AbstractUnits = constructorof(U)(dimension(u), remove_offset(u.todims))
is_scalar(u::AbstractUnits) = is_scalar(todims(u))
udynamic(u::Units) = Units(udynamic(dimension(u)), todims(u), usymbol(u))
ustatic(u::Units) = Units(ustatic(dimension(u)), todims(u), usymbol(u))
dimtype(::Type{U}) where {D,U<:Units{D}} = dimtype(D)
dimvaltype(::Type{U}) where {D,U<:Units{D}} = dimvaltype(D)
dimval(u::Units) = dimval(dimension(u))
unit(x::NumUnion) = Units(NoDims(), NoTransform(), DEFAULT_USYMBOL)

#=================================================================================================
Quantity types
The basic type is Quantity, which belongs to <:Any (hence it has no real hierarchy)
=================================================================================================#
"""
    struct Quantity{T<:Number, U<:AbstractUnitLike} <: Number
        value :: T
        unit  :: U
    end

Numeric quantity type (that sutypes to number) with fields `value` and `unit`
"""
struct Quantity{T<:Number, U<:AbstractUnitLike} <: Number
    value :: T
    unit  :: U
end

"""
    struct FlexQuant{T<:Any, U<:AbstractUnitLike}
        value :: T
        unit  :: U
    end

Generic quantity type (that can hold any value) with fields `value` and `unit`
"""
struct FlexQuant{T<:Any, U<:AbstractUnitLike}
    value :: T
    unit  :: U
end

"""
    QuantUnion{T<:Any,U<:AbstractUnitLike}

Convenience union that allows Number and Non-Number types to be considered together
"""
const QuantUnion{T,U} = Union{Quantity{T,U}, FlexQuant{T,U}}

Quantity{T}(x, u::AbstractUnitLike) where T = Quantity{T, typeof(u)}(x, u)
#Quantity(x::T, u::StaticUnits{D}) where {T,D} = Quantity(u.todims(x), StaticDims{D}())
#Quantity{T}(x, u::StaticUnits{D}) where {T,D} = Quantity{T}(convert(T, u.todims(x)), StaticDims{D}())
Quantity{T}(q::QuantUnion) where T = Quantity{T}(ustrip(q), unit(q))
Quantity{T,U}(q::QuantUnion) where {T,U} = Quantity{T,U}(ustrip(q), unit(q))
Quantity{T,StaticDims{d}}(q::QuantUnion) where {T,d} = Quantity{T,StaticDims{d}}(ustrip(d, q), StaticDims{d}())
Quantity{T,StaticDims{d}}(x::Number) where {T,d} = Quantity{T, StaticDims{d}}(convert(T,x), StaticDims{d}())
Quantity{T,D}(q::QuantUnion) where {T,D<:AbstractDimensions} = Quantity{T,D}(dstrip(q), dimension(q))
Quantity{T,D}(x::Number) where {T,D<:AbstractDimensions} = Quantity{T,D}(x, D())

FlexQuant{T}(x, u::AbstractUnitLike) where T = FlexQuant{T, typeof(u)}(x, u)
#FlexQuant(x::T, u::StaticUnits{D}) where {T,D} = FlexQuant(u.todims(x), StaticDims{D}())
#FlexQuant{T}(x, u::StaticUnits{D}) where {T,D} = FlexQuant{T}(convert(T, u.todims(x)), StaticDims{D}())
FlexQuant{T}(q::QuantUnion) where T = FlexQuant{T}(ustrip(q), unit(q))
FlexQuant{T,U}(q::QuantUnion) where {T,U} = FlexQuant{T,U}(ustrip(q), unit(q))

ustrip(q::QuantUnion) = q.value
unit(q::QuantUnion) = q.unit
todims(q::QuantUnion) = todims(q.unit)
dimension(q::QuantUnion) = dimension(unit(q))
unittype(::Type{<:QuantUnion{T,U}}) where {T,U} = U
dimtype(::Type{<:QuantUnion{T,U}}) where {T,U} = dimtype(U)
dimvaltype(::Type{<:QuantUnion{T,U}}) where {T,U} = dimvaltype(U)
udynamic(q::QuantUnion) = Quantity(ustrip(q), udynamic(unit(q)))
valtype(q::Any) = valtype(typeof(q))
valtype(::Type{T}) where T = T
valtype(::Type{<:QuantUnion{T}}) where T = T 


#Translation between AffineTansform and numeric quantities
AffineTransform(scale::Real, offset::Quantity) = AffineTransform(scale=scale, offset=dstrip(offset))

#Generic promotion functions of quantities based on element type 
"""
    quantity(x, u::AbstractUnitLike)

Constructs a quantity with arguments `x` and `u`, and selects the appropriate type to use based on the arguments
(if x is a number, `Quantity<:Number` is used, otherwise `FlexQuant<:Any` is used)
"""
quantity(x, u::AbstractUnitLike) = quant_type(x)(x, u)

"""
    quant_type(::Type{T})

Returns the appropriate no-parameter quantity TypeUnion associated with type `T`
(if `T <: Number`, `Quantity` is returned, otherwise `FlexQuant` is returned)
"""
quant_type(::Type{<:Number}) = Quantity 
quant_type(::Type{<:Any}) = FlexQuant
quant_type(x) = quant_type(typeof(x))



"""
    constructorof(::Type{T}) where T = Base.typename(T).wrapper

Return the constructor of a type `T{PS...}` by default it only returns `T` (i.e. removes type parameters)
This function can be overloaded if custom behaviour is needed
"""
constructorof(::Type{T}) where T = Base.typename(T).wrapper
constructorof(::Type{<:Dimensions}) = Dimensions
constructorof(::Type{<:Units}) = Units
constructorof(::Type{<:Quantity}) = Quantity
constructorof(::Type{<:FlexQuant}) = FlexQuant


#=============================================================================================
Errors and assertion functions
=============================================================================================#

"""
    DimensionError{D} <: Exception

Error thrown when an operation is dimensionally invalid given the arguments
"""
struct DimensionError{T} <: Exception
    items :: T
end
DimensionError(arg1, arg2, args...) = DimensionError((arg1, arg2, args...))
Base.showerror(io::IO, e::DimensionError{<:Tuple}) = print(io, "DimensionError: ", e.items, " have incompatible dimensions")
Base.showerror(io::IO, e::DimensionError{<:QuantUnion}) = print(io, "DimensionError: ", e.items, " is not dimensionless")
Base.showerror(io::IO, e::DimensionError{<:AbstractUnitLike}) = print(io, "DimensionError: ", e.items, " is not dimensionless")

"""
    ConversionError{U, U0} <: Exception 

Error thrown when trying to convert u0 to u, offers hint on how to make u compatible
"""
struct ConversionError{U<:AbstractUnitLike, U0<:AbstractUnitLike} <: Exception
    u::U
    u0::U0
end

function Base.showerror(io::IO, e::ConversionError{<:AbstractUnitLike, <:AbstractUnitLike})
    io_tmp = IOBuffer()
    pretty_str(x) = (ushow(io_tmp, x, pretty=true); String(take!(io_tmp)))

    uΔ = dimension(e.u0)/dimension(e.u)
    return print(io, "ConversionError: Cannot convert unit '", pretty_str(e.u0), "' to target unit '", pretty_str(e.u), "'. Consider multiplying '", pretty_str(e.u), "' by '", pretty_str(uΔ), "' or similar.")
end

"""
    NotScalarError{D} <: Exception

Error thrown for non-scalar units when the operation is only valid for scalar units
"""
struct NotScalarError{D} <: Exception
    dim::D
    NotScalarError(dim) = new{typeof(dim)}(dim)
end
Base.showerror(io::IO, e::NotScalarError) = print(io, "NotScalarError: ", e.dim, " cannot be treated as scalar, operation only valid for scalar units")


"""
    NotDimensionError{D} <: Exception

Error thrown for non-dimensional units (scaled or affine) when the operation is only valid for dimensional units
"""
struct NotDimensionError{D} <: Exception
    dim::D
    NotDimensionError(dim) = new{typeof(dim)}(dim)
end
Base.showerror(io::IO, e::NotDimensionError) = print(io, "NotDimensionError: ", e.dim, " cannot be treated as dimension, operation only valid for dimension units")


assert_scalar(u::AbstractDimLike) = u
assert_scalar(u::AbstractUnits) = is_scalar(u) ? u : throw(NotScalarError(u))
scalar_dimension(u::AbstractUnitLike) = dimension(assert_scalar(u))

assert_dimension(u::AbstractDimLike) =  u
assert_dimension(u::AbstractUnits) = is_dimension(u) ? u : throw(NotDimensionError(u))

assert_dimensionless(u::AbstractUnitLike) = isdimensionless(u) ? u : throw(DimensionError(u))
assert_dimensionless(q::QuantUnion) = isdimensionless(unit(q)) ? q : throw(DimensionError(q))
dimensionless(u::AbstractUnitLike) = dimension(assert_dimensionless(u))
dimensionless(q::QuantUnion) = ustrip(assert_dimensionless(ubase(q)))
dimensionless(n) = n

isdimensionless(u::AbstractUnitLike) = isdimensionless(dimension(u))
isdimensionless(::Type{StaticDims{d}}) where d = isdimensionless(d)
isdimensionless(::Type{Units{StaticDims{d},T}}) where {d,T} = isdimensionless(d)
isdimensionless(d::AbstractDimLike)  = iszero(d) || isunknown(d)
isdimensionless(d::NoDims) = true
Base.iszero(u::D) where D<:AbstractDimensions = (u == D(0))
Base.iszero(u::StaticDims{d}) where d = iszero(d)