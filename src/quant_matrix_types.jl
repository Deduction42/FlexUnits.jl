using LinearAlgebra
using StaticArrays

import StaticArrays.StaticLUMatrix

const UnitOrDims{D} = Union{D, AbstractUnits{D}} where D<:AbstractDimensions
const ScalarOrVec{T} = Union{T, AbstractVector{T}} where T

abstract type AbstractUnitMap{U<:UnitOrDims{<:AbstractDimensions}} <: AbstractMatrix{U} end

"""
    QuantMatrixUnits{U<:UnitOrDims, A<:AbstractMatrix{<:QuantUnion{<:Any,U}}}

Wraps a quantity matrix A so that getting its index returns the units of the element
"""
struct QuantMatrixUnits{U<:UnitOrDims, A<:AbstractMatrix{<:QuantUnion{<:Any,U}}} <: AbstractMatrix{U}
    quantmat :: A
    QuantMatrixUnits(a::A) where {U,A<:AbstractMatrix{QuantUnion{<:Any,U}}} = new{U,A}(a)
end
Base.IndexStyle(::Type{QuantMatrixUnits{D,A}}) where{D,A} = IndexStyle(A)
Base.getindex(m::QuantMatrixUnits, args...) = broadcast(unit, getindex(m.quantmat, args...))

"""
    QuantMatrixDims{D<:AbstractDimensions, A<:AbstractMatrix{<:QuantUnion{<:Any,<:UnitOrDims{D}}}}

Wraps a quantity matrix A so that getting its index returns the dimensiosn of the element
"""
struct QuantMatrixDims{D<:AbstractDimensions, A<:AbstractMatrix{<:QuantUnion{<:Any,<:UnitOrDims{D}}}} <: AbstractMatrix{D}
    quantmat :: A
    QuantMatrixDims(a::A) where {D, A<:AbstractMatrix{<:QuantUnion{<:Any,<:UnitOrDims{D}}}} = new{D,A}(a)
end
Base.IndexStyle(::Type{QuantMatrixDims{D,A}}) where{D,A} = IndexStyle(A)
Base.getindex(m::QuantMatrixDims, args...) = broadcast(dimension, getindex(m.quantmat, args...))

Base.size(m::Union{<:QuantMatrixUnits,<:QuantMatrixDims}) = size(m.quantmat)

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

function UnitMap(md::AbstractMatrix{<:AbstractDimensions})
    u_out = md[:,begin]
    u_in = u_out[begin]./md[begin,:]

    #Check for dimensional consistency
    for jj in axes(md,2), ii in axes(md,1)
        (u_out[ii]/md[ii,jj] == u_in[jj]) || error("Unit inconsistency around index $([ii,jj]) of original matrix, expected dimension '$(u_out[ii]*inv(u_in[jj]))', found unit '$(unit(md[ii,jj]))'")
    end
    return UnitMap(u_out=u_out, u_in=u_in)
end

UnitMap(mq::AbstractMatrix{<:QuantUnion}) = UnitMap(QuantMatrixDims(mq))



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
struct LinmapQuant{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:UnitMaps{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    units :: U
end

A linear mapping of quantities. A special kind of matrix that is intended to be used for multiplying vectors of quantities;
such matrices must be dimensionally consistent and can be represented by a UnitMap. These constraints lead to much faster 
unit inference and a smaller memory footprint (M+N instead of M*N for units).
"""
struct LinmapQuant{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:AbstractUnitMap{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    units :: U
end

function LinmapQuant(::Type{U}, m::AbstractMatrix{<:Quantity}) where U <: AbstractUnitMap
    values = ustrip_base.(m)
    units  = U(m)
    return LinmapQuant(values, units)
end

LinmapQuant(m::AbstractMatrix{<:Quantity}) = LinmapQuant(UnitMap, m)
ustrip(lq::LinmapQuant) = lq.values
unit(lq::LinmapQuant) = lq.units

Base.getindex(q::LinmapQuant, ii::Integer, jj::Integer) = q.values[ii,jj] * q.units[ii,jj]
Base.size(q::LinmapQuant) = size(q.values)
Base.inv(q::LinmapQuant) = LinmapQuant(inv(q.values), inv(q.units))

#======================================================================================================================
Linear Mapping Factorizations
======================================================================================================================#
"""
struct FactorQuant{T, D<:AbstractDimensions, F<:Factorization{T}, U<:AbstractUnitMap{D}}
    factor :: F
    units  :: U 
end

A factored linear mapping. A subclass of Factorizations with a unit mapping attached. Calling getproperty
is re-routed to the original factor, with the appropriate units calcualted from the mapping.
"""
struct FactorQuant{F, D<:AbstractDimensions, U<:AbstractUnitMap{D}}
    factor :: F
    units  :: U 
end
ustrip(fq::FactorQuant) = getfield(fq, :factor)
unit(fq::FactorQuant) = getfield(fq, :units)
Base.inv(fq::FactorQuant) = LinmapQuant(inv(ustrip(fq)), inv(unit(fq)))
LinearAlgebra.inv!(fq::FactorQuant) = LinmapQuant(inv!(ustrip(fq)), inv(unit(fq)))

# LU Factorization ===================================================================================
LinearAlgebra.lu(mq::LinmapQuant; kwargs...) = FactorQuant(lu(ustrip(mq); kwargs...), unit(mq))

#May need to iterate over more subtypes of AbstractMatrix
LinearAlgebra.lu(mq::AbstractMatrix{<:Quantity}, args...; kwargs...) = lu(LinmapQuant(UnitMap, mq), args...; kwargs...)
StaticArrays.lu(mq::StaticLUMatrix{N,M,<:Quantity}; kwargs...) where {N,M} = lu(LinmapQuant(UnitMap, mq); kwargs...)

function Base.getproperty(fq::FactorQuant{<:LU, D}, fn::Symbol) where D
    F = ustrip(fq)

    if fn === :L
        u = unit(fq)
        return LinmapQuant(F.L, UnitMap(u_in=uinput(u).^0, u_out=uoutput(u)[invperm(F.p)]))
    elseif fn === :U 
        u = unit(fq)
        return LinmapQuant(F.L, UnitMap(u_in=uinput(u), u_out=uoutput(u).^0))
    elseif fn === :p 
        return F.p
    elseif fn === :P
        P0 = F.P 
        T  = eltype(P0)
        Pu = zeros(Quantity{T,D}, size(P0))
        for (ii, x) in enumerate(P0)
            isone(x) && setindex!(Pu, x, ii)
        end
        return Pu
    else
        getfield(fq, fn)
    end
end

#======================================================================================================================
Nonlinear mapping
======================================================================================================================#
"""
struct FunctionQuant{F, U<:AbstractUnitMap}
    func  :: F
    units :: U
end

A generic mapping with units. Useful for applying units to unitless functions that assume units for inputs/outputs.
"""
struct FunctionQuant{F, U<:AbstractUnitMap}
    func  :: F
    units :: U
end

function (qmap::FunctionQuant)(x)
    fmap = qmap.func
    umap = qmap.units
    xraw = _strictmap(ustrip, uinput(umap), x)
    return _strictmap(*, fmap(xraw), uoutput(umap))
end

function _strictmap(f, args...)
    allequal(map(length, args)) || throw(ArgumentError("All arguments must be of equal length"))
    return map(f, args...)
end




