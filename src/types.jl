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
    AbstractDimLike

A broad class representing anything that can be interpreted as a dimension
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


const UnitOrDims{D} = Union{D, AbstractUnits{D}} where D<:AbstractDimLike
const ScalarOrVec{T} = Union{T, AbstractVector{T}} where T

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
uscale(u::AbstractUnitLike) = uscale(tobase(u))
uoffset(u::AbstractUnitLike) = uoffset(tobase(u))
dimtype(x::Any) = dimtype(typeof(x))
dimtype(::Type{T}) where T = error("dimtype not yet implemented for type $(T)")
dimtype(::Type{D}) where D<:AbstractDimLike = D
dimtype(::Type{U}) where {D, U<:AbstractUnits{D}} = dimtype(D)
dimvaltype(x::Any) = dimvaltype(typeof(x)) 
dimvaltype(::Type{T}) where T = dimtype(T)
dimpowtype(::Type{D}) where {P, D<:AbstractDimensions{P}} = P
Base.getindex(d::AbstractDimensions, k::Symbol) = getproperty(d, k)
dimpowtype(::Type{U}) where {U<:AbstractUnitLike} = dimpowtype(dimtype(U))
is_dimension(u::AbstractUnitLike) = is_identity(tobase(u))
is_scalar(u::AbstractUnitLike) = is_scalar(tobase(u))
udynamic(u::AbstractDimensions) = u

#=======================================================================================
Transforms
=======================================================================================#
"""
NoTransform object, the default transform returned by tobase(x::AbstractDimensionLike). Calling it results in 
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
Base.inv(t::NoTransform, x) = x
uscale(t::NoTransform)  = 1
uoffset(t::NoTransform) = 0
tobase(u::AbstractDimLike) = NoTransform()
is_identity(t::NoTransform) = true
is_scalar(t::NoTransform) = true


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
(t::AffineTransform)(t0::AbstractUnitTransform) = t ∘ t0

(t::AffineTransform)(x) = x*t.scale + t.offset
Base.inv(t::AffineTransform, x) = (x - t.offset)/t.scale
(t::AffineTransform)(x::AbstractArray) = t.(x)
Base.inv(t::AffineTransform, x::AbstractArray) = inv.(t, x)
(t::AffineTransform)(x::Tuple) = map(t, x)
Base.inv(t::AffineTransform, x::Tuple) = map(xi->inv(t, xi), x)

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
    @kwdef struct ExpAffTransform{T<:Real} <: AbstractUnitTransform
        scale  :: T = 1.0
        offset :: T = 0.0
    end

A type representing an affine transfomration formula that applies "exp" at the end;
this is useful for logarithmic units

# Constructors
    ExpAffTransform(scale::Real, offset::Real)
    ExpAffTransform(; scale, offset)
"""

@kwdef struct ExpAffTransform{T<:Real} <: AbstractUnitTransform
    scale  :: T = 1.0
    offset :: T = 0.0
end
ExpAffTransform(scale::T1, offset::T2) where {T1,T2} = ExpAffTransform{promote_type(T1,T2)}(scale, offset)
(t::ExpAffTransform)(t0::AbstractUnitTransform) = t ∘ t0

uscale(t::ExpAffTransform) = t.scale 
uoffset(t::ExpAffTransform) = t.offset

(t::ExpAffTransform)(x) = exp(x*t.scale + t.offset)
Base.inv(t::ExpAffTransform, x) = (log(x) - t.offset)/t.scale
(t::ExpAffTransform)(x::AbstractArray) = t.(x)
Base.inv(t::ExpAffTransform, x::AbstractArray) = inv.(t, x)
(t::ExpAffTransform)(x::Tuple) = map(t, x)
Base.inv(t::ExpAffTransform, x::Tuple) = map(xi->inv(t, xi), x)

Base.inv(t::ExpAffTransform) = Base.Fix1(inv, t)

#=======================================================================================
Dimensions
=======================================================================================#
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
Dimensions(args::Real...) = Dimensions{FixRat32}(args...)

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
(::Type{D})(d::StaticDims) where {D<:AbstractDimensions} = D(dimval(d))

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
Units
=======================================================================================#
"""
    @kwdef struct Units{D<:AbstractDimLike, T<:AbstractUnitTransform} <: AbstractUnits{D, T}
        dims   :: D
        tobase :: T
        symbol :: Symbol = DEFAULT_USYMBOL
    end

A dynamic unit object that contains dimensions (dims) and its conversion formula to said dimensions (tobase). The conversion 
formula determines what kind of unit is referred to. An AffineTransform implies affine units, a NoTransform implies dimensions.
Units with dynamic dimensions can generated through the `@ud_str` macro, while units with static dimensions can be generated
throiugh the `@u_str` macro.

# Constructors
    Units{D}(units, tobase::T, symbol=DEFAULT_USYMBOL) where {D,T<:AbstractUnitTransform} 
    Units(dims::D, tobase::AbstractUnitTransform=NoTransform(), symbol=DEFAULT_USYMBOL) where D<:AbstractDimensions 
    Units(units::D, tobase::AbstractUnitTransform, symbol=DEFAULT_USYMBOL) where D<:AbstractUnits

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
    tobase :: T
    symbol :: Symbol = DEFAULT_USYMBOL
end
Units{D}(dims, tobase::T, symbol=DEFAULT_USYMBOL) where {D,T<:AbstractUnitTransform} = Units{D,T}(dims, tobase, symbol)
Units(dims::D, tobase::AbstractUnitTransform=NoTransform(), symbol=DEFAULT_USYMBOL) where D<:AbstractDimLike = Units(dims, tobase, symbol)
Units(units::U, tobase::AbstractUnitTransform, symbol=DEFAULT_USYMBOL) where U<:AbstractUnits = Units(dimension(assert_dimension(units)), tobase, symbol)

tobase(u::Units) = u.tobase
dimension(u::Units) = u.dims 
usymbol(u::Units) = u.symbol
remove_offset(u::U) where U<:AbstractUnits = constructorof(U)(dimension(u), remove_offset(u.tobase))
is_scalar(u::AbstractUnits) = is_scalar(tobase(u))
udynamic(u::Units) = Units(udynamic(dimension(u)), tobase(u), usymbol(u))
ustatic(u::Units) = Units(ustatic(dimension(u)), tobase(u), usymbol(u))
dimtype(::Type{U}) where {D,U<:Units{D}} = dimtype(D)
dimvaltype(::Type{U}) where {D,U<:Units{D}} = dimvaltype(D)
dimval(u::Units) = dimval(dimension(u))
dimval(::Type{<:Units{StaticDims{D}}}) where D = D
unit(x::NumUnion) = Units(NoDims(), NoTransform(), DEFAULT_USYMBOL)

#=================================================================================================
Quantities
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

Quantity{T,U}(q::QuantUnion) where {T,U} = Quantity{T,U}(ustrip(q), unit(q))
Quantity{T,D}(q::QuantUnion) where {T,D<:StaticDims} = Quantity{T,D}(ustrip(D(), q), D())
Quantity{T,D}(q::QuantUnion) where {T,D<:AbstractDimensions} = Quantity{T,D}(dstrip(q), dimension(q))

Quantity{T,U}(x::NumUnion) where {T,U} = Quantity{T,U}(x*dimtype(U)())
Quantity{T,D}(x::NumUnion) where {T,D<:StaticDims} = Quantity{T,D}(convert(T,x), assert_dimensionless(D()))
Quantity{T,D}(x::NumUnion) where {T,D<:AbstractDimensions} = Quantity{T,D}(x, D())

Quantity{T}(x, u::AbstractUnitLike) where T = Quantity{T, typeof(u)}(x, u)
Quantity{T}(q::QuantUnion) where T = Quantity{T}(ustrip(q), unit(q))

FlexQuant{T,U}(q::QuantUnion) where {T,U} = FlexQuant{T,U}(ustrip(q), unit(q))
FlexQuant{T,D}(q::QuantUnion) where {T,D<:StaticDims} = FlexQuant{T,D}(ustrip(D(), q), D())
FlexQuant{T,D}(q::QuantUnion) where {T,D<:AbstractDimensions} = Quantity{T,D}(dstrip(q), dimension(q))

FlexQuant{T}(x, u::AbstractUnitLike) where T = FlexQuant{T, typeof(u)}(x, u)
FlexQuant{T}(q::QuantUnion) where T = FlexQuant{T}(ustrip(q), unit(q))

"""
    ubase(q::QuantUnion)

Converts quantity `q` to its raw dimensional equivalent (such as SI units)
"""
function ubase(q::QuantUnion{<:Any,<:AbstractUnitLike})
    u  = unit(q)
    ft = tobase(u)
    return quantity(ft(ustrip(q)), dimension(u))
end 
ubase(q::QuantUnion{<:Any,<:AbstractDimLike}) = q

ustrip(q::QuantUnion) = q.value
unit(q::QuantUnion) = q.unit
tobase(q::QuantUnion) = tobase(q.unit)
dstrip(q::QuantUnion) = tobase(q)(ustrip(q))
ustrip_base(q::QuantUnion) = dstrip(q)
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


#=================================================================================================
Logarithmic quantities (only supports numeric types)
=================================================================================================#
"""
    struct LogQuant{T<:Number, U<:AbstractUnitLike} <: Number 
        value :: T 
        unit  :: U 
    end

Numeric quantity type (that subtypes to number) representing the lograrithm of a Quantity.
This changes arithmeic operations (+-*/) according to logarithmic algebra. 
"""
struct LogQuant{T<:Number, U<:AbstractUnitLike} <: Number 
    value :: T 
    unit  :: U 
end
LogQuant{T}(x, u::AbstractUnitLike) where T = LogQuant{T, typeof(u)}(x, u)
LogQuant{T}(q::QuantUnion) where T = LogQuant{T}(log(dstrip(q)), dimension(q))
LogQuant{T,U}(q::QuantUnion) where {T,U} = LogQuant{T,U}(log(dstrip(q)), dimension(q))

Quantity{T,U}(lq::LogQuant) where {T,U} = Quantity{T,U}(ubase(lq))
FlexQuant{T,U}(lq::LogQuant) where {T,U} = FlexQuant{T,U}(ubase(lq))

const LogQuantUnion{T,U} = Union{Quantity{T,U}, FlexQuant{T,U}, LogQuant{T,U}}

logquant(x::T, u::AbstractUnitLike) where T = LogQuant{T}(x, u)
logquant(x::T, u::Units{<:AbstractDimLike, <:ExpAffTransform}) where {T} = logquant(x, Units(dimension(u), log(tobase(u)), usymbol(u)))
logquant(q::Quantity) = LogQuant(log(dstrip(q)), dimension(q))
logquant(lq::LogQuant) = lq 
quantity(q::LogQuant) = ubase(q)
linquant(q::Quantity) = q
linquant(q::LogQuant) = ubase(q)

"""
    logubase(lq::LogQuant)

Converts log-quantity `lq` to its natural logarithmic scale (Nepers)
"""
function logubase(lq::LogQuant{<:Any,<:AbstractUnitLike})
    u  = unit(lq)
    ft = tobase(u)
    return logquant(ft(ustrip(lq)), dimension(u))
end 
logubase(lq::LogQuant{<:Any,<:AbstractDimLike}) = lq
logubase(q::Quantity) = log(ubase(q))

function ubase(q::LogQuant{<:Any,<:AbstractUnitLike})
    u  = unit(q)
    ft = exp(tobase(u))
    return quantity(ft(ustrip(q)), dimension(u))
end 

ustrip(q::LogQuant) = q.value
unit(q::LogQuant) = q.unit
tobase(q::LogQuant) = tobase(q.unit)
dstrip(q::LogQuant) = exp(tobase(q)(ustrip(q)))
ustrip_base(q::LogQuant) = dstrip(q)
dimension(q::LogQuant) = dimension(unit(q))
unittype(::Type{<:LogQuant{T,U}}) where {T,U} = U
dimtype(::Type{<:LogQuant{T,U}}) where {T,U} = dimtype(U)
dimvaltype(::Type{<:LogQuant{T,U}}) where {T,U} = dimvaltype(U)
udynamic(q::LogQuant) = Quantity(ustrip(q), udynamic(unit(q)))
valtype(::Type{<:LogQuant{T}}) where T = T 

#=================================================================================================
Logarithmic scale objects used to scale LogQuant
=================================================================================================#
"""
    struct LogScale{T<:Real}
        scale :: T
        base :: T
        symbol :: Symbol
    end

A callable object used to apply a scale to a logarithmic unit. For example,
```julia
dB = LogScale(scale=0.1, base=10, symbol=:dB)
```
After defining it, this object can also be called on a LogQuant to change the logarithmic scale
```julia
julia> dB(log(10u"kPa"))
40.0 dB(kg/(m s²))
```
It can also be used to create logarithmic units that you can convert to
julia> 10u"kPa" |> dB(u"Pa")
39.99999999999999 dB(Pa)
"""
@kwdef struct LogScale{T<:Real}
    scale :: T
    base :: T
    symbol :: Symbol
end
LogScale(scale::T1, base::T2, symbol) where {T1,T2} = LogScale{promote_type(T1, T2)}(scale, base, symbol)

function (s::LogScale)(reference::Union{AbstractUnitLike, Quantity})
    return Units(
        dims = dimension(reference),
        tobase = ExpAffTransform(
            scale = s.scale*log(s.base), 
            offset = log(dstrip(1*reference)),
        ),
        symbol = (s.symbol==DEFAULT_USYMBOL) ? s.symbol : Symbol(string(s.symbol)*"($(reference))")
    )
end

(s::LogScale)(q::LogQuant) = uconvert(s(dimension(q)), q)
(s::LogScale)() = s(NoDims())

dB = LogScale(scale=0.1, base=10, symbol=:dB)
Np = LogScale(scale=1, base=exp(1), symbol=:Np)


#=============================================================================================
Constructor utilities
=============================================================================================#
"""
    constructorof(::Type{T}) where T = Base.typename(T).wrapper

Return the constructor of a type `T{PS...}` by default it only returns `T` (i.e. removes type parameters)
This function can be overloaded if custom behaviour is needed
"""
constructorof(::Type{T}) where T = Base.typename(T).wrapper
constructorof(::Type{<:Dimensions}) = Dimensions
constructorof(::Type{<:Units}) = Units
constructorof(::Type{<:Quantity}) = Quantity


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
    return print(io, "ConversionError: Cannot convert unit '",pretty_str(e.u0),"' to target unit '",pretty_str(e.u),"' due to a dimension mismatch of '",pretty_str(uΔ),"'")
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
Base.showerror(io::IO, e::NotDimensionError) = print(io, "NotDimensionError: ", e.dim," cannot be treated as dimension, operation only valid for dimension units")

"""
    LogLinearError{F} <: Exception

Error thrown for function {F} when applied to a mixture of linear and logaritmic units
"""
struct LogLinearError{F, Q1, Q2} <: Exception
    f :: F
    q1 :: Q1 
    q2 :: Q2 
end
Base.showerror(io::IO, e::LogLinearError{F,Q1,Q2}) where {F,Q1,Q2} = print(io, 
    "LogLinearError: Cannot apply function '",e.f,"' to argument types ",Q1," and ",Q2,
    "'. Perhaps you meant to convert one of the arguments to linear units using 'ubase(x::LoqQuant)'"
)

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