using LinearAlgebra
using StaticArrays

import StaticArrays.StaticLUMatrix

const UnitOrDims{D} = Union{D, AbstractUnits{D}} where D<:AbstractDimensions
const ScalarOrVec{T} = Union{T, AbstractVector{T}} where T

abstract type AbstractUnitMap{U<:UnitOrDims{<:AbstractDimensions}} end
abstract type AbstractDimsMap{D<:AbstractDimensions} <: AbstractMatrix{D} end


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
unit(m::AbstractMatirx) = QuantMatrixUnits(m)

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
dimension(m::AbstractMatrix) = QuantMatrixDims(m)

Base.size(m::Union{<:QuantMatrixUnits,<:QuantMatrixDims}) = size(m.quantmat)

"""
struct UnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}, TO<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_in  :: TI
    u_out :: TO
end

Used to represent a unit transformation from input units 'u_in' to outpout units 'u_out'. Often applied
to nonlinear functions.
"""
@kwdef struct UnitMap{U<:UnitOrDims, TI<:ScalarOrVec{U}, TO<:ScalarOrVec{U}} <: AbstractUnitMap{U}
    u_in  :: TI
    u_out :: TO
end


"""
struct DimsMap{D<:AbstractDimLike, TI<:AbstractVector{D}, TO<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_in  :: TI
    u_out :: TO
end

Used to represent a unit transformation from input dimensions 'u_in' to outpout dimensions 'u_out'.
This is like a unit map but focuses on dimensions, simplifying linear algebra (it subtypes to Matrix)
"""
@kwdef struct DimsMap{D<:AbstractDimensions, TI<:AbstractVector{D}, TO<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_in  :: TI
    u_out :: TO
end

Base.getindex(m::DimsMap, ii::Integer, jj::Integer) = m.u_out[ii]/m.u_in[jj]
Base.size(m::DimsMap) = (length(m.u_out), length(m.u_in))
Base.inv(m::DimsMap) = DimsMap(u_out=m.u_in, u_in=m.u_out)
Base.adjoint(m::DimsMap) = DimsMap(u_out=inv.(m.u_in), u_in=inv.(m.u_out))
Base.transpose(m::DimsMap) = adjoint(m)
uoutput(m::DimsMap) = m.u_out
uinput(m::DimsMap) = m.u_in

function DimsMap(md::AbstractMatrix{<:AbstractDimensions})
    u_out = md[:,begin]
    u_in = u_out[begin]./md[begin,:]

    #Check for dimensional consistency
    for jj in axes(md,2), ii in axes(md,1)
        (u_out[ii]/md[ii,jj] == u_in[jj]) || error("Unit inconsistency around index $([ii,jj]) of original matrix, expected dimension '$(u_out[ii]*inv(u_in[jj]))', found unit '$(unit(md[ii,jj]))'")
    end
    return DimsMap(u_out=u_out, u_in=u_in)
end

DimsMap(mq::AbstractMatrix{<:QuantUnion}) = DimsMap(QuantMatrixDims(mq))

function canonical!(u::DimsMap)
    ui1 = u.u_in[1]
    unew = if isdimensionless(ui1)
        u
    elseif ArrayInterface.can_setindex(u.u_in) && ArrayInterface.can_setindex(u.u_out)
        u.u_in  .= u.u_in ./ ui1
        u.u_out .= u.uout ./ ui1
        u
    else
        DimsMap(
            u_in = u.in./ui1, 
            u_out = u_out./ui1
        )
    end
    return unew
end



"""
struct RepDimsMap{D<:Abstractdimensions, TI<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_scale :: D
    u_in :: TI
end

Used to represent a special kind of dimensional transformation (a repeatable transformation) of dimensions 'u_in'.
If "U" is repeatable, "U*U*...*U*x" is a valid operation and the units of "U*x" are similar to the units of "x".
If 'u_scale' is dimensionless, the unit transformation is idempotent (same output units as input). 
This structure enables certain kinds of operations such as matrix powers (whose unit transform must be repeatable).
Idempotence enables even more transformations like matrix exponentials.
"""
@kwdef struct RepDimsMap{D<:AbstractDimensions, TI<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_scale :: D
    u_in :: TI
end

Base.getindex(m::RepDimsMap, ii::Integer, jj::Integer) = m.u_scale*m.u_in[ii]/m.u_in[jj]
Base.size(m::RepDimsMap) = (length(m.u_in), length(m.u_in))
Base.inv(m::RepDimsMap)  = RepDimsMap(u_scale=inv(m.u_scale), u_in=m.u_in)
Base.adjoint(m::RepDimsMap) = RepDimsMap(u_scale=m.u_scale, u_in=inv.(m.u_in))
uoutput(m::RepDimsMap) = map(Base.Fix1(*, m.u_in), m.u_scale)
uinput(m::RepDimsMap) = m.u_in

function RepDimsMap(md::DimsMap)
    #Matrix must be square
    sz = size(md)
    sz[1] == sz[2] || throw(DimensionMismatch("Repeatable Unit Mapping must be square: dimensions are $(sz)"))

    #Calculate the uniform scale
    u_scale = md.u_out[begin]/md.u_in[begin]

    #Verify that uniform scale is consistent 
    for (u_out, u_in) in zip(md.u_out, md.u_in)
        u_out/u_in == u_scale || error("Cannot convert to Repeatable Unit Mapping: $(md.u_out) and $(md.u_in) must share a common factor")
    end

    return RepDimsMap(u_scale=u_scale, u_in=md.u_in)
end

function RepDimsMap(mq::AbstractMatrix{<:Union{QuantUnion,AbstractUnitLike}})
    return RepDimsMap(DimsMap(mq))
end

DimsMap(md::RepDimsMap) = DimsMap(u_out=md.u_in.*md.u_scale, u_in=md.u_in)

function canonical!(u::RepDimsMap)
    ui1 = u.u_in[1]
    u_in = if isdimensionless(ui1)
        u.u_in
    elseif ArrayInterface.can_setindex(u.u_in)
        u.u_in .= u.u_in ./ ui1
        u.u_in
    else
        u.in./ui1
    end
    return RepDimsMap(u_in = u_in, u_scale = u.u_scale/ui1)
end


"""
struct SymUnitMap{D<:AbstractDimensions, TI<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_scale :: D
    u_in :: TI
end

Used to represent a special kind of dimensional transformation (a symmetric transformation) of dimensions 'u_in'.
If "U" is symmetric, then "x'U*x" is a valid operation and the dimensions of "U*x" are similar to the inverse of the
dimensions of "x". This structure enables certain kinds of operations reserved for symmetric matrices.
"""
@kwdef struct SymUnitMap{D<:AbstractDimensions, TI<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_scale :: D
    u_in :: TI
end

Base.getindex(m::SymUnitMap, ii::Integer, jj::Integer) = m.u_scale*m.u_in[ii]*m.u_in[jj]
Base.size(m::SymUnitMap) = (length(m.u_in), length(m.u_in))
Base.inv(m::SymUnitMap)  = SymUnitMap(u_scale=inv(m.u_scale), u_in=m.u_in)
Base.adjoint(m::SymUnitMap) = SymUnitMap(u_scale=m.u_scale, u_in=inv.(m.u_in))
uoutput(m::SymUnitMap) = map(Base.Fix1(*, m.u_uscale), m.u_in)
uinput(m::SymUnitMap) = m.u_in

function SymUnitMap(md::DimsMap)
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
    return SymUnitMap(DimsMap(mq))
end

DimsMap(md::SymUnitMap) = DimsMap(u_out=md.u_in.*md.u_scale, u_in=md.u_in)

function canonical!(u::SymUnitMap)
    ui1  = u.u_in[1]
    u_in = if isdimensionless(ui1)
        u.u_in
    elseif ArrayInterface.can_setindex(u.u_in)
        u.u_in .= u.u_in ./ ui1
        u.u_in
    else
        u.in./ui1
    end
    return SymUnitMap(u_in = u_in, u_scale = u.u_scale*ui1)
end

"""
struct LinmapQuant{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:UnitMaps{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    units :: U
end

A linear mapping of quantities. A special kind of matrix that is intended to be used for multiplying vectors of quantities;
such matrices must be dimensionally consistent and can be represented by a UnitMap. These constraints lead to much faster 
unit inference and a smaller memory footprint (M+N instead of M*N for units).
"""
struct LinmapQuant{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:AbstractDimsMap{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    dims :: U
end

function LinmapQuant(m::AbstractMatrix{T}, u::UnitMap) where T 
    todims(u::AbstractUnits, n) = u.todims(n)
    new_m = todims.(m, u.u_out./u.u_in')
    new_u = DimsMap(u_in = dimension.(u.u_in), u_out = dimension.(u.u_out))
    return LinmapQuant(new_m, canonical!(new_u))
end

function LinmapQuant(::Type{U}, m::AbstractMatrix{<:Quantity}) where U <: AbstractDimsMap
    values = dstrip.(m)
    units  = U(dimension(m))
    return LinmapQuant(values, units)
end
LinmapQuant(m::AbstractMatrix{<:Quantity}) = LinmapQuant(DimsMap, m)

ustrip(lq::LinmapQuant) = lq.values
unit(lq::LinmapQuant) = lq.dims
dimension(lq::LinmapQuant) = lq.dims
ubase(lq::LinmapQuant) = lq

Base.getindex(q::LinmapQuant, ii::Integer, jj::Integer) = q.values[ii,jj] * q.dims[ii,jj]
Base.size(q::LinmapQuant) = size(q.values)
Base.inv(q::LinmapQuant) = LinmapQuant(inv(q.values), inv(q.dims))

#======================================================================================================================
Linear Mapping Factorizations
======================================================================================================================#
"""
struct FactorQuant{T, D<:AbstractDimensions, F<:Factorization{T}, U<:AbstractUnitMap{D}}
    factor :: F
    dims  :: U 
end

A factored linear mapping. A subclass of Factorizations with a unit mapping attached. Calling getproperty
is re-routed to the original factor, with the appropriate units calcualted from the mapping.
"""
struct FactorQuant{F, D<:AbstractDimensions, U<:AbstractDimsMap{D}}
    factor :: F
    dims  :: U 
end
ustrip(fq::FactorQuant) = getfield(fq, :factor)
unit(fq::FactorQuant) = getfield(fq, :dims)
dimension(fq::FactorQuant) = getfield(fq, :dims)
Base.inv(fq::FactorQuant) = LinmapQuant(inv(ustrip(fq)), inv(unit(fq)))
LinearAlgebra.inv!(fq::FactorQuant) = LinmapQuant(inv!(ustrip(fq)), inv(unit(fq)))

# LU Factorization ===================================================================================
LinearAlgebra.lu(mq::LinmapQuant; kwargs...) = FactorQuant(lu(ustrip(mq); kwargs...), unit(mq))

#May need to iterate over more subtypes of AbstractMatrix
LinearAlgebra.lu(mq::AbstractMatrix{<:Quantity}, args...; kwargs...) = lu(LinmapQuant(DimsMap, mq), args...; kwargs...)
StaticArrays.lu(mq::StaticLUMatrix{N,M,<:Quantity}; kwargs...) where {N,M} = lu(LinmapQuant(DimsMap, mq); kwargs...)

function Base.getproperty(fq::FactorQuant{<:LU, D}, fn::Symbol) where D
    F = ustrip(fq)

    if fn === :L
        u = unit(fq)
        return LinmapQuant(F.L, DimsMap(u_in=uinput(u).^0, u_out=uoutput(u)[invperm(F.p)]))
    elseif fn === :U 
        u = unit(fq)
        return LinmapQuant(F.L, DimsMap(u_in=uinput(u), u_out=uoutput(u).^0))
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
Eigenvalue decomposition will have a slightly different approach

1.  Firstly, the eigenvalues will always be dimensionless 
2.  Secondly, the eigenvector matrix will have a dimensionless output, but dimensionful input 
    -   Can be performed on both repeatable and symmetric matrices
        -   this will require "issymmetric" and "isrepeatable" functions
    -   Symmetric matrices are recovered using the transpose of the DimsMap
    -   Repeatable matrices are recovered using the inverse of the DimsMap
======================================================================================================================#

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
struct FunctionQuant{F, U<:AbstractDimsMap}
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


#======================================================================================================================
Special cases
======================================================================================================================#
#UniformScaling with dynamic dimensions should produce unknown dimension on off-diagonals (consistent with other behaviour)
Base.getindex(J::UniformScaling{T}, i::Integer, j::Integer) where T<:Quantity = ifelse(i==j, J.λ, zero(T))

#=
for op in (:+, :-)
    @eval function Base.$op(m::AbstractMatrix{<:QuantUnion{<:Any, <:AbstractDimensions}}, s::UniformScaling{Quantity{<:Any, <:StaticDims}})
        return $op(m, UniformScaling(udynamic(s[1,1])))
    end

    @eval function Base.$op(s::UniformScaling{Quantity{<:Any, <:StaticDims}}, m::AbstractMatrix{<:QuantUnion{<:Any, <:AbstractDimensions}})
        return $op(m, UniformScaling(udynamic(s[1,1])))
    end
end
=#

#=
Recently found errors where adding the following produces an error:

Matrix{Quantity{<:Any, <:AbstractDimensions}} + UniformScaling{Quantity{<:Any, <:StaticDims}}

If the matrix isn't uniform units, the addiiton fails, because all off-diagonals have known units if the quantity is static 
We may want to promote the uniform scaling in this case to be dynamic, so that off-diagonals have unknown units
=#