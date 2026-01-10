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

function UnitMap(mq::AbstractMatrix{<:Union{QuantUnion,AbstractUnitLike}})
    u_out = dimension.(mq[:,begin])
    u_in = u_out[begin]./dimension.(mq[begin,:])

    #Check for dimensional consistency
    for jj in axes(mq,2), ii in axes(mq,1)
        (u_out[ii]/dimension(mq[ii,jj]) == u_in[jj]) || error("Unit inconsistency around index $([ii,jj]) of original matrix, expected dimension '$(u_out[ii]*inv(u_in[jj]))', found unit '$(unit(mq[ii,jj]))'")
    end

    return UnitMap(u_out=u_out, u_in=u_in)
end

"""
struct RepUnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_scale :: U
    u_in :: TI
end

Used to represent a special kind of unit transformation (a repeatable transformation) of units 'u_in'.
If "U" is repeatable, "U*U*...*U*x" is a valid operation and the units of "U*x" are similar to the units of "x".
If 'u_scale' is dimensionless, the unit transformation is idempotent (same output units as input). 
This structure enables certain kinds of operations such as matrix powers (whose unit transform must be repeatable).
Idempotence enables even more transformations like matrix exponentials.
"""
@kwdef struct RepUnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_scale :: U
    u_in :: TI
end

Base.getindex(m::RepUnitMap, ii::Integer, jj::Integer) = m.u_scale*m.u_in[ii]/m.u_in[jj]
Base.size(m::RepUnitMap) = (length(m.u_in), length(m.u_in))
Base.inv(m::RepUnitMap)  = RepUnitMap(u_scale=inv(m.u_scale), u_in=m.u_in)
Base.adjoint(m::RepUnitMap) = RepUnitMap(u_scale=m.u_scale, u_in=inv.(m.u_in))
uoutput(m::RepUnitMap) = map(Base.Fix1(*, m.u_in), m.u_scale)
uinput(m::RepUnitMap) = m.u_in

function RepUnitMap(md::UnitMap)
    #Matrix must be square
    sz = size(md)
    sz[1] == sz[2] || throw(DimensionMismatch("Repeatable Unit Mapping must be square: dimensions are $(sz)"))

    #Calculate the uniform scale
    u_scale = md.u_out[begin]/md.u_in[begin]

    #Verify that uniform scale is consistent 
    for (u_out, u_in) in zip(md.u_out, md.u_in)
        u_out/u_in == u_scale || error("Cannot convert to Repeatable Unit Mapping: $(md.u_out) and $(md.u_in) must share a common factor")
    end

    return RepUnitMap(u_scale=u_scale, u_in=md.u_in)
end

function RepUnitMap(mq::AbstractMatrix{<:Union{QuantUnion,AbstractUnitLike}})
    return RepUnitMap(UnitMap(mq))
end

UnitMap(md::RepUnitMap) = UnitMap(u_out=md.u_in.*md.u_scale, u_in=md.u_in)



"""
struct SymUnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_scale :: U
    u_in :: TI
end

Used to represent a special kind of unit transformation (a symmetric transformation) of units 'u_in'.
If "U" is symmetric, then "x'U*x" is a valid operation and the units of "U*x" are similar to the inverse units of "x".
This unit structure enables certain kinds of operations reserved for symmetric matrices.
"""
@kwdef struct SymUnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_scale :: U
    u_in :: TI
end

Base.getindex(m::SymUnitMap, ii::Integer, jj::Integer) = m.u_scale*m.u_in[ii]*m.u_in[jj]
Base.size(m::SymUnitMap) = (length(m.u_in), length(m.u_in))
Base.inv(m::SymUnitMap)  = SymUnitMap(u_scale=inv(m.u_scale), u_in=m.u_in)
Base.adjoint(m::SymUnitMap) = SymUnitMap(u_scale=m.u_scale, u_in=inv.(m.u_in))
uoutput(m::SymUnitMap) = map(Base.Fix1(*, m.u_uscale), m.u_in)
uinput(m::SymUnitMap) = m.u_in

function SymUnitMap(md::UnitMap)
    #Matrix must be square
    sz = size(md)
    sz[1] == sz[2] || throw(DimensionMismatch("Symmetric Unit Mapping must be square: dimensions are $(sz)"))

    #Calculate the uniform scale
    u_scale = md.u_out[begin]*md.u_in[begin]

    #Verify that the uniform scale is consistent
    for (u_out, u_in) in zip(md.u_out, md.u_in)
        u_out*u_in == u_scale || error("Cannot convert to Symmetric Unit Mapping: $(md.u_in) and $(md.u_out) must be similar inverses")
    end

    return SymUnitMap(u_scale=u_scale, u_in=inv.(md.u_in))
end

function SymUnitMap(mq::AbstractMatrix{<:Union{QuantUnion,AbstractUnitLike}})
    return SymUnitMap(UnitMap(mq))
end

UnitMap(md::SymUnitMap) = UnitMap(u_out=md.u_in.*md.u_scale, u_in=md.u_in)



"""
struct QuantLinMap{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:UnitMaps{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    units :: U
end

A linear mapping of quantities. A special kind of matrix that is intended to be used for multiplying vectors of quantities;
such matrices must be dimensionally consistent and can be represented by a UnitMap. These constraints lead to much faster 
unit inference and a smaller memory footprint (M+N instead of M*N for units).
"""
struct QuantLinMap{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:AbstractUnitMap{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    units :: U
end

function QuantLinMap(::Type{U}, m::AbstractMatrix{<:Quantity}) where U <: AbstractUnitMap
    values = ustrip_base.(m)
    units  = U(dimension.(m))
    return QuantLinMap(values, units)
end

QuantLinMap(m::AbstractMatrix{<:Quantity}) = QuantLinMap(UnitMap, m)

Base.getindex(q::QuantLinMap, ii::Integer, jj::Integer) = q.values[ii,jj] * q.units[ii,jj]
Base.size(q::QuantLinMap) = size(q.values)
Base.inv(q::QuantLinMap) = QuantLinMap(inv(q.values), inv(q.units))


#======================================================================================================================
Nonlinear mapping
======================================================================================================================#
"""
struct QuantMapping{F, U<:AbstractUnitMap}
    func  :: F
    units :: U
end

A generic mapping with units. Useful for applying units to unitless functions that assume units for inputs.
"""
struct QuantMapping{F, U<:AbstractUnitMap}
    func  :: F
    units :: U
end

function (qmap::QuantMapping)(x)
    fmap = qmap.func
    umap = qmap.units
    xraw = _strictmap(ustrip, uinput(umap), x)
    return _strictmap(*, fmap(xraw), uoutput(umap))
end

function _strictmap(f, args...)
    allequal(map(length, args)) || throw(ArgumentError("All arguments must be of equal length"))
    return map(f, args...)
end

#======================================================================================================================
Linear algebra relationships with "Matrix" 
The outer type specializes first, so something like
    Base.inv(q::AbstractMatrix{<:Quantity}) = inv(QuantTransform(DimTransform, q))
Will be skipped in the case of Matrix{<:Quantity} (it will use inv(m::Matrix) in Base instead)
Because of this, such code will have to specify all desired concrete matrix types
======================================================================================================================#
Base.inv(q::Matrix{<:Quantity}) = inv(QuantLinMap(UnitMap, q))


#======================================================================================================================
Shortcut multiplication strategies
*(U::DimTransform, M::AbstractMatrix) = U.u_out * (U.u_in'*M)
*(M::AbstractMatrix, U::DimTransform) = (M*U.u_out) * U.u_in'
*(U1::DimTransform, U2::DimTransform) = U1.u_out.*dot(U1.u_in, U2.u_out).*U2.u_in'
      = DimTransform(u_out=U1.u_out.*dot(U1.u_in, U2.u_out), u_in=U2.u_in) 
*(U1::ScaleDimTransform, U2::ScaleDimTransform) = DimTransform(u_scale=U1.u_scale*U2.u_scale, u_out=U1.u_out) iff U1.u_out == U2.u_out
======================================================================================================================#
using Test
import .UnitRegistry.@u_str
using StaticArrays
import Random
using Statistics

#Nonlinear map
@kwdef struct PumpInput{T} <: FieldVector{2,T}
    current :: T 
    voltage :: T
end

@kwdef struct PumpOutput{T} <: FieldVector{3,T}
    power :: T 
    pressure :: T
    flow :: T 
end

function pumpfunc(x::PumpInput)
    p = x.current*x.voltage*0.9   
    return PumpOutput(power = p, pressure = sqrt(p), flow = sqrt(p))
end
pumpfunc(x::AbstractVector) = pumpfunc(PumpInput(x))


@testset "Linear Mapping Basics" begin
    Random.seed!(1234)

    #Quck tests 
    u1 = [u"lbf*ft", u"kW", u"rpm"]
    u2 = [u"kg/s", u"m^3/hr", u"kW"]

    xm = randn(3,3)
    qM = QuantLinMap(UnitMap, xm.*u2./u1')

    #Matrix inversion
    x = randn(3).*u1
    y = qM*x
    @test all(x .≈ inv(qM)*y)

    #Matrix transpose
    @test all(Matrix(y') .≈ Matrix(x'*qM'))

    #Square matrices
    Σ = cov(randn(20,3)*rand(3,3))
    x = randn(3).*u2

    #Symmetric matrix
    rS = Σ.*inv.(u2).*inv.(u2)'
    qS = QuantLinMap(SymUnitMap, rS)
    @test all(qS .≈ ubase.(rS))
    @test x'*(rS)*x ≈ x'*qS*x

    #Repeatable matrix
    rR = Σ.*u2.*inv.(u2)'
    qR = QuantLinMap(RepUnitMap, rR)
    @test all(qR .≈ ubase.(rR))
    @test all(rR^2*x .≈ qR^2*x)

    #Nonlinear mapping
    pumpunits = UnitMap(PumpInput(current=u"A", voltage=u"V"), PumpOutput(power=u"W", pressure=u"Pa", flow=u"m^3/s"))
    upumpfunc = QuantMapping(pumpfunc, pumpunits)
    qinput = PumpInput(current=500*u"mA", voltage=6u"V")
    @test all(upumpfunc(qinput) .≈ pumpfunc(ustrip.(uinput(pumpunits), qinput)).*uoutput(pumpunits))

end