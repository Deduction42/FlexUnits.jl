#Preamble (delete when finished)
include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")

const UnitOrDims{D} = Union{D, AbstractUnits{D}} where D<:AbstractDimensions
const ScalarOrVec{T} = Union{T, AbstractVector{T}} where T

abstract type AbstractUnitMap{U<:UnitOrDims{<:AbstractDimensions}} <: AbstractMatrix{U} end

"""
struct UnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}, TO<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_in  :: TI
    u_out :: TO
end

Used to represent a unit transformation from input units 'u_in' to outpout units 'u_out'.
This is often applied to matrices or functions of vectors.
"""
@kwdef struct UnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}, TO<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_in  :: TI
    u_out :: TO
end

Base.getindex(m::UnitMap, ii::Integer, jj::Integer) = m.u_out[ii]/m.u_in[jj]
Base.size(m::UnitMap) = (length(m.u_out), length(m.u_in))
Base.inv(m::UnitMap) = UnitMap(u_out=m.u_in, u_in=m.u_out)
Base.adjoint(m::UnitMap) = UnitMap(u_out=inv.(m.u_in), u_in=inv.(m.u_out))
Base.transpose(m::UnitMap) = adjoint(m)
uoutput(m::UnitMap) = m.u_out
uinput(m::UnitMap) = m.u_in

function UnitMap(mq::AbstractMatrix{<:Union{AbstractQuantity,AbstractUnitLike}})
    u_out = dimension.(mq[:,begin])
    u_in = u_out[begin]./dimension.(mq[begin,:])

    #Check for dimensional consistency
    for jj in axes(mq,2), ii in axes(mq,1)
        (u_out[ii]/dimension(mq[ii,jj]) == u_in[jj]) || error("Unit inconsistency around index $([ii,jj]) of original matrix, expected dimension '$(u_out[ii]*inv(u_in[jj]))', found unit '$(unit(mq[ii,jj]))'")
    end

    return UnitMap(u_out=u_out, u_in=u_in)
end

"""
struct RUnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_scale :: U
    u_in :: TI
end

Used to represent a special kind of unit transformation (a recursive or repeatable transformation) of units 'u_in'.
If "U" is recursive, "U*U*...*U*x" is a valid operation; under such circumstances, the units of "U*x" are the same as "x" 
ith a potential uniform unit scale, "u_scale". If 'u_scale' is dimensionless, the unit transformation is idempotent (same output 
units as input). This structure enables certain kinds of operations such as matrix powers (whose unit transform must be recursive).
Idempotence enables even more transformations like matrix exponentials.
"""
@kwdef struct RUnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_scale :: U
    u_in :: TI
end

Base.getindex(m::RUnitMap, ii::Integer, jj::Integer) = m.u_scale*m.u_in[ii]/m.u_in[jj]
Base.size(m::RUnitMap) = (length(m.u_in), length(m.u_in))
Base.inv(m::RUnitMap)  = RUnitMap(u_scale=inv(m.u_scale), u_in=m.u_in)
Base.adjoint(m::RUnitMap) = RUnitMap(u_scale=m.u_scale, u_in=inv.(m.u_in))
uoutput(m::RUnitMap) = map(Base.Fix1(*, m.u_in), m.u_scale)
uinput(m::RUnitMap) = m.u_in

function RUnitMap(md::UnitMap)
    #Matrix must be square
    sz = size(md)
    sz[1] == sz[2] || throw(DimensionMismatch("Recursive Unit Mapping must be square: dimensions are $(sz)"))

    #Calculate the uniform scale
    u_scale = md.u_out[begin]/md.u_in[begin]

    #Verify that uniform scale is consistent 
    for (u_out, u_in) in zip(md.u_out, md.u_in)
        u_out/u_in == u_scale || error("Cannot convert to Recursive Unit Mapping: $(md.u_out) and $(md.u_in) must share a common factor")
    end

    return RUnitMap(u_scale=u_scale, u_in=md.u_in./u_scale)
end

function RUnitMap(mq::AbstractMatrix{<:Union{AbstractQuantity,AbstractUnitLike}})
    return RUnitMap(UnitMap(mq))
end

UnitMap(md::RUnitMap) = UnitMap(u_out=md.u_in.*md.u_scale, u_in=inv.(md.u_in))

const UnitMaps{U,V} = Union{UnitMap{U,V}, RUnitMap{U,V}}

"""
struct LinMapQuant{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:UnitMaps{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    units :: U
end

A linear mapping quantity. A special kind of matrix that is intended to be used for multiplying vectors of quantities;
such matrices must be dimensionally consistent and can be represented by a UnitMap. These constraints lead to much faster 
unit inference and a smaller memory footprint (M+N instead of M*N for units).
"""
struct LinMapQuant{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:UnitMaps{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    units :: U
end

function LinMapQuant(::Type{U}, m::AbstractMatrix{<:Quantity}) where U <: AbstractUnitMap
    values = ustrip_base.(m)
    units  = U(dimension.(m))
    return LinMapQuant(values, units)
end

LinMapQuant(m::AbstractMatrix{<:Quantity}) = LinMapQuant(UnitMap, m)

Base.getindex(q::LinMapQuant, ii::Integer, jj::Integer) = q.values[ii,jj] * q.units[ii,jj]
Base.size(q::LinMapQuant) = size(q.values)
Base.inv(q::LinMapQuant) = LinMapQuant(inv(q.values), inv(q.units))


#======================================================================================================================
Linear algebra relationships with "Matrix" 
The outer type specializes first, so something like
    Base.inv(q::AbstractMatrix{<:Quantity}) = inv(QuantTransform(DimTransform, q))
Will be skipped in the case of Matrix{<:Quantity} (it will use inv(m::Matrix) in Base instead)
Because of this, such code will have to specify all desired concrete matrix types
======================================================================================================================#
Base.inv(q::Matrix{<:Quantity}) = inv(LinMapQuant(UnitMap, q))


#======================================================================================================================
Shortcut multiplication strategies
*(U::DimTransform, M::AbstractMatrix) = U.u_out * (U.u_in'*M)
*(M::AbstractMatrix, U::DimTransform) = (M*U.u_out) * U.u_in'
*(U1::DimTransform, U2::DimTransform) = U1.u_out.*dot(U1.u_in, U2.u_out).*U2.u_in'
      = DimTransform(u_out=U1.u_out.*dot(U1.u_in, U2.u_out), u_in=U2.u_in) 
*(U1::ScaleDimTransform, U2::ScaleDimTransform) = DimTransform(u_scale=U1.u_scale*U2.u_scale, u_out=U1.u_out) iff U1.u_out == U2.u_out
======================================================================================================================#
import .UnitRegistry.@u_str

#Quck tests 
u1 = [u"lbf*ft", u"kW", u"rpm"]
u2 = [u"kg/s", u"m^3/hr", u"kW"]

xm = randn(3,3)
qm = LinMapQuant(UnitMap, xm.*u2./u1')

x = randn(3).*u1
y = qm*x
inv(qm)*y

