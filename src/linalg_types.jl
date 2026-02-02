using LinearAlgebra
using StaticArrays
using SparseArrays

import ArrayInterface
import StaticArrays.StaticLUMatrix
import SparseArrays.AbstractSparseMatrixCSC
import SparseArrays.AbstractCompressedVector

dimtype(::Type{<:AbstractUnitMap{U}}) where U = dimtype(U)
dimvaltype(::Type{<:AbstractUnitMap{U}}) where U = dimvaltype(U)

# AbstractDimsMap is similar to AbstractArray, but doesn't subtype to it; it subtypes to AbstractUnitMap which is not really an array
# Moreover subtyping to AbstractArray produces ambiguity problems, it's easier to just redefine the API methods
abstract type AbstractDimsMap{D<:AbstractDimensions} <: AbstractUnitMap{D} end

Base.eltype(::Type{<:AbstractDimsMap{D}}) where D = D
Base.IndexStyle(::Type{<:AbstractDimsMap}) = IndexCartesian()
Base.getindex(m::AbstractDimsMap, ind::CartesianIndex{2}) = m[ind[1],ind[2]]
Base.getindex(m::AbstractDimsMap, ind::Integer) = m[CartesianIndices(m)[ind]]
#Base.iterate(m::AbstractDimsMap, i=1) = (@inline; (i - 1)%UInt < length(m)%UInt ? (m[i], i + 1) : nothing)
Base.CartesianIndices(m::AbstractDimsMap) = CartesianIndices(axes(m))
Base.length(m::AbstractDimsMap) = prod(size(m))
Base.collect(m::AbstractDimsMap) = uoutput(m).*inv.(uinput(m)')
Base.axes(m::AbstractDimsMap, d::Integer) = d <= 2 ? axes(m)[d] : OneTo(1)
assert_symmetric(m::AbstractDimsMap)  = issymmetric(m) ? m : throw(ArgumentError("Input must be symmetric. Recieved $(m)"))
assert_repeatable(m::AbstractDimsMap) = isrepeatable(m) ? m : throw(ArgumentError("Input must support repeatable multiplication. Received $(m)"))
assert_idempotent(m::AbstractDimsMap) = isidempotent(m) ? m : throw(ArgumentError("Input must be idempotent. Received $(m)"))
LinearAlgebra.issymmetric(m::AbstractDimsMap) = all(d[1]==inv(d[2]) for d in strictzip(uinput(m), uoutput(m)))
isrepeatable(m::AbstractDimsMap) = all(d[1]==d[2] for d in strictzip(uinput(m), uoutput(m)))
isidempotent(m::AbstractDimsMap) = isdimensionless(m.u_fac) && isrepeatable(m)


"""
    QuantArrayDims{D<:AbstractDimensions, A<:AbstractArray} <: AbstractArray{D}

Wraps a quantity matrix A so that getting its index returns the dimensiosn of the element
"""
struct QuantArrayDims{D<:AbstractDimensions, N, A<:AbstractArray} <: AbstractArray{D,N}
    array :: A
    function QuantArrayDims(a::AbstractArray{<:Any,N}) where N
        D = dimvaltype(eltype(a))
        return new{D,N,typeof(a)}(a)
    end
end
Base.IndexStyle(::Type{<:QuantArrayDims{D,N,A}}) where{D,N,A} = IndexStyle(A)
Base.getindex(m::QuantArrayDims, args...) = broadcast(dimension, getindex(m.array, args...))
Base.size(m::QuantArrayDims) = size(m.array)
Base.axes(m::QuantArrayDims) = axes(m.array)
ArrayInterface.can_setindex(::Type{QuantArrayDims}) = false
dimension(m::AbstractArray) = QuantArrayDims(m)
dimension(m::SArray) = dimension.(m)


"""
    QuantArrayVals{D<:AbstractDimensions, A<:AbstractArray} <: AbstractArray{D}

Wraps a quantity matrix A so that getting its index returns the value of hte element
"""
struct QuantArrayVals{T, N, A<:AbstractArray} <: AbstractArray{T,N}
    array :: A
    function QuantArrayVals(a::AbstractArray{<:Any,N}) where N
        T = valtype(eltype(a))
        return new{T,N,typeof(a)}(a)
    end
end
Base.IndexStyle(::Type{<:QuantArrayVals{D,N,A}}) where{D,N,A} = IndexStyle(A)
Base.getindex(m::QuantArrayVals, args...) = broadcast(dstrip, getindex(m.array, args...))
Base.size(m::QuantArrayVals) = size(m.array)
Base.axes(m::QuantArrayVals) = axes(m.array)
ArrayInterface.can_setindex(::Type{QuantArrayVals}) = false
dstrip(m::AbstractArray) = QuantArrayVals(m)
dstrip(m::SArray) = dstrip.(m)

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
uoutput(m::UnitMap) = m.u_out
uinput(m::UnitMap) = m.u_in


"""
struct DimsMap{D<:AbstractDimLike, TI<:AbstractVector{D}, TO<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_fac :: D
    u_in  :: TI
    u_out :: TO
end

Used to represent a unit transformation from input dimensions 'u_in' to outpout dimensions 'u_out'.
This is like a unit map but focuses on dimensions and has matrix-like behaviour since dimensions 
support linear algebra, but generic units may not (affine units, logarithmic units etc).

WARNING: DimsMap constructor expects u_in and u_out to be scaled so that the first element is zero.
In order to prevent excessive allocations, u_in and u_out will be mutated in-place if possible to fit these
requirements. Ensure these requirements are met, supply immutable arguments or copies if mutating is undesirable.
"""
@kwdef struct DimsMap{D<:AbstractDimensions, TI<:AbstractVector{D}, TO<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_fac :: D
    u_in  :: TI
    u_out :: TO
    function DimsMap{D, TI, TO}(u_fac, u_in_raw, u_out_raw) where {D<:AbstractDimensions, TI<:AbstractVector{D}, TO<:AbstractVector{D}}
        (u_fac, u_in)  = canonical_input!(u_fac, u_in_raw)
        (u_fac, u_out) = canonical_output!(u_fac, u_out_raw)
        return new{D, TI, TO}(u_fac, u_in, u_out)
    end
    function DimsMap(u_fac::D, u_in::TI, u_out::TO) where {D<:AbstractDimensions, TI<:AbstractVector{D}, TO<:AbstractVector{D}}
        return DimsMap{D, TI, TO}(u_fac, u_in, u_out)
    end
end

Base.axes(m::DimsMap) = (axes(m.u_out)[1], axes(m.u_in)[1])
Base.size(m::DimsMap) = (length(m.u_out), length(m.u_in))
Base.inv(m::DimsMap) = DimsMap(u_fac=inv(m.u_fac), u_out=m.u_in, u_in=m.u_out)
uoutput(m::DimsMap) = m.u_out
uinput(m::DimsMap) = m.u_in
ufactor(m::DimsMap) = m.u_fac
Base.getindex(m::DimsMap, ii::Integer, jj::Integer) = m.u_out[ii]/m.u_in[jj]*m.u_fac
Base.getindex(m::DimsMap, ii::Integer, vj::Any) = (m.u_fac*m.u_out[ii]) ./ m.u_in[vj]
Base.getindex(m::DimsMap, vi::Any, jj::Integer) = (m.u_fac/m.u_in[jj]) .* m.u_out[vi]
Base.getindex(m::DimsMap, vi::Any, vj::Any) = DimsMap(u_fac=m.u_fac, u_out=m.u_out[vi], u_in=m.u_in[vj])

function DimsMap(md::AbstractMatrix{<:AbstractDimensions})
    u_fac = md[begin,begin]
    u_out = md[:,begin]./u_fac
    u_in  = u_fac./md[begin,:]

    #Check for dimensional consistency
    for jj in axes(md,2), ii in axes(md,1)
        md_ij = u_out[ii]/u_in[jj]*u_fac
        (md_ij == md[ii, jj]) || error("Unit inconsistency around index $([ii, jj]) of original matrix, expected dimension '$(md_ij))', found dimension '$(md[ii, jj])'")
    end
    return DimsMap(u_fac=u_fac, u_out=u_out, u_in=u_in)
end

DimsMap(mq::AbstractMatrix{<:QuantUnion}) = DimsMap(QuantArrayDims(mq))
DimsMap(d::AbstractDimsMap) = d


"""
    AdjointDmap{D, M<:AbstractUnitMap{D}} <: AbstractUnitMap{D}

Wraps a dimension map as an adjoint/transpose (and avoids subtyping to AbstractArray)
"""
struct AdjointDmap{D, M<:AbstractDimsMap{D}} <: AbstractDimsMap{D}
    parent :: M 
end
Base.IndexStyle(::Type{AdjointDmap{D,M}}) where {D,M} = IndexStype(M)
Base.transpose(m::AbstractDimsMap) = AdjointDmap(m)
Base.transpose(m::AdjointDmap) = m.parent
Base.adjoint(m::AbstractDimsMap) = AdjointDmap(m)
Base.adjoint(m::AdjointDmap) = m.parent 
Base.axes(m::AdjointDmap) = reverse(axes(m.parent))
Base.inv(m::AdjointDmap) = adjoint(inv(m.parent))
uoutput(m::AdjointDmap) = inv.(uinput(m.parent))
uinput(m::AdjointDmap) = inv.(uoutput(m.parent))
ufactor(m::AdjointDmap) = ufactor(m.parent)
LinearAlgebra.issymmetric(m::AdjointDmap) = issymmetric(m.parent)
isrepeatable(m::AdjointDmap) = isrepeatable(m.parent)
isidempotent(m::AdjointDmap) = isidempotent(m.parent)
Base.getindex(m::AdjointDmap, ind1::Integer, ind2::Integer) = getindex(m.parent, ind2, ind1)
Base.getindex(m::AdjointDmap, ind1::Any, ind2::Any) = AdjointDmap(getindex(m.parent, ind2, ind1))


"""
struct LinmapQuant{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:UnitMaps{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    units :: U
end

A linear mapping of quantities. A special kind of matrix that is intended to be used for multiplying vectors of quantities;
such matrices must be dimensionally consistent and can be represented by a UnitMap. These constraints lead to much faster 
unit inference and a smaller memory footprint (O(M+N) instead of O(M*N*N2) in the case of multiplication).
"""
struct LinmapQuant{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:AbstractDimsMap{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    dims :: U
end

function LinmapQuant(m::AbstractMatrix{T}, u::UnitMap) where T 
    todims(u::AbstractUnits, n) = u.todims(n)
    todims(u::AbstractDimensions, n) = n
    new_m = todims.(u.u_out./u.u_in', m) 
    new_u = DimsMap(u_fac = zero(dimvaltype(u)), u_in = dimension.(u.u_in), u_out = dimension.(u.u_out))
    return LinmapQuant(new_m, new_u)
end

LinmapQuant(m::QuantArrayVals, d::QuantArrayDims) = LinmapQuant(dstrip.(m.array), DimsMap(d))
LinmapQuant(m::AbstractMatrix) = LinmapQuant(dstrip.(m), DimsMap(dimension(m)))
LinmapQuant(m::LinmapQuant) = m

ustrip(lq::LinmapQuant) = lq.values
dstrip(lq::LinmapQuant) = lq.values
unit(lq::LinmapQuant) = lq.dims
dimension(lq::LinmapQuant) = lq.dims
ubase(lq::LinmapQuant) = lq

Base.IndexStyle(::Type{<:LinmapQuant}) = IndexCartesian()
Base.size(q::LinmapQuant) = size(q.values)
Base.inv(q::LinmapQuant) = LinmapQuant(inv(q.values), inv(q.dims))
Base.transpose(q::LinmapQuant) = LinmapQuant(transpose(q.values), transpose(q.dims))
Base.adjoint(q::LinmapQuant) = LinmapQuant(adjoint(q.values), adjoint(q.dims))
Base.getindex(q::LinmapQuant, ii::Integer, jj::Integer) = q.values[ii,jj] * q.dims[ii,jj]
Base.getindex(q::LinmapQuant, ii::Integer, vj::Any) = VectorQuant(q.values[ii,vj], q.dims[ii,vj])
Base.getindex(q::LinmapQuant, vi::Any, jj::Integer) = VectorQuant(q.values[vi,jj], q.dims[vi,jj])
Base.getindex(q::LinmapQuant, vi::Any, vj::Any) = LinmapQuant(q.values[vi,vj], q.dims[vi,vj])

"""
struct VectorQuant{T, D<:AbstractDimensions, V<:AbstractVector{T}, U<:AbstractVector{D}} <: AbstractVector{Quantity{T,D}}
    values :: V
    units :: U
end

A vector of quantities. A special kind of vector that separates dimensions from values for easier dimension manipulation when used
with LinmapQuant
"""
struct VectorQuant{T, D<:AbstractDimensions, V<:AbstractVector{T}, U<:AbstractVector{D}} <: AbstractVector{Quantity{T,D}}
    values :: V
    dims :: U
end

function VectorQuant(v::AbstractVector{T}, u::AbstractVector{<:AbstractUnits}) where T 
    todims(ux::AbstractUnits, x) = ux.todims(x)
    new_m = todims.(v, u)
    new_u = dimension.(u)
    return VectorQuant(new_m, new_u)
end

VectorQuant(m::QuantArrayVals, d::QuantArrayDims) = VectorQuant(dstrip.(m.array), dimension.(d.array))
VectorQuant(v::AbstractVector) = VectorQuant(dstrip.(v), dimension.(v))
VectorQuant(v::VectorQuant) = v

ustrip(lq::VectorQuant) = lq.values
dstrip(lq::VectorQuant) = lq.values
unit(lq::VectorQuant) = lq.dims
dimension(lq::VectorQuant) = lq.dims
ubase(lq::VectorQuant) = lq

Base.IndexStyle(::Type{<:VectorQuant{<:Any, <:Any, V}}) where {V} = IndexStyle(V)
Base.getindex(q::VectorQuant, ii::Integer) = q.values[ii] * q.dims[ii]
Base.getindex(q::VectorQuant, vi::Any) = VectorQuant(q.values[vi], q.dims[vi])
Base.size(q::VectorQuant) = size(q.values)

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
dstrip(fq::FactorQuant) = getfield(fq, :factor)
unit(fq::FactorQuant) = getfield(fq, :dims)
dimension(fq::FactorQuant) = getfield(fq, :dims)
Base.inv(fq::FactorQuant) = LinmapQuant(inv(ustrip(fq)), inv(unit(fq)))
LinearAlgebra.inv!(fq::FactorQuant) = LinmapQuant(inv!(ustrip(fq)), inv(unit(fq)))
Base.transpose(fq::FactorQuant) = FactorQuant(transpose(fq.factor), transpose(fq.dims))
Base.adjoint(fq::FactorQuant) = FactorQuant(adjoint(fq.factor), adjoint(fq.dims))

# LU Factorization ===================================================================================
LinearAlgebra.lu(mq::LinmapQuant; kwargs...) = FactorQuant(lu(ustrip(mq); kwargs...), unit(mq))
LinearAlgebra.lu(mq::LinmapQuant, ::Val{true}; kwargs...) = FactorQuant(lu(ustrip(mq), Val(true); kwargs...), unit(mq))
LinearAlgebra.lu(mq::LinmapQuant, ::Val{false}; kwargs...) = FactorQuant(lu(ustrip(mq), Val(false); kwargs...), unit(mq))

#May need to iterate over more subtypes of AbstractMatrix
StaticArrays.lu(mq::StaticLUMatrix{N,M,<:Quantity}; kwargs...) where {N,M} = lu(LinmapQuant(mq); kwargs...)
LinearAlgebra.lu(mq::AbstractMatrix{<:Quantity}; kwargs...) = lu(LinmapQuant(mq); kwargs...)
LinearAlgebra.lu(mq::AbstractMatrix{<:Quantity}, ::Val{true}; kwargs...) = lu(LinmapQuant(mq), Val(true); kwargs...)
LinearAlgebra.lu(mq::AbstractMatrix{<:Quantity}, ::Val{false}; kwargs...) = lu(LinmapQuant(mq), Val(false); kwargs...)

function Base.getproperty(fq::FactorQuant{<:Union{LU, StaticArrays.LU}, D}, fn::Symbol) where D
    F = ustrip(fq)

    if fn === :L
        u = unit(fq)
        return LinmapQuant(F.L, DimsMap(u_fac=ufactor(u)^0, u_in=uinput(u).^0, u_out=uoutput(u)[F.p]))
    elseif fn === :U 
        u = unit(fq)
        return LinmapQuant(F.U, DimsMap(u_fac=ufactor(u), u_in=uinput(u), u_out=uoutput(u).^0))
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
struct FunctionQuant{F, U<:AbstractUnitMap}
    func  :: F
    units :: U
end

function (qmap::FunctionQuant)(x)
    fmap = qmap.func
    umap = qmap.units
    xraw = strictmap(ustrip, uinput(umap), x)
    return strictmap(*, fmap(xraw), uoutput(umap))
end

function strictmap(f, args...)
    allequal(map(length, args)) || throw(ArgumentError("All arguments must be of equal length"))
    return map(f, args...)
end

function strictzip(args...)
    allequal(map(length, args)) || throw(ArgumentError("All arguments must be of equal length"))
    return zip(args...)
end

#======================================================================================================================
Special cases
======================================================================================================================#
#UniformScaling with dynamic dimensions should produce unknown dimension on off-diagonals (consistent with other behaviour)
Base.getindex(J::UniformScaling{T}, i::Integer, j::Integer) where T<:Quantity = ifelse(i==j, J.λ, zero(T))


#======================================================================================================================
Utility functions
======================================================================================================================#
function canonical_input!(u_fac::D, u_in::V) where {D<:AbstractDimensions, V<:AbstractVector{D}}
    u0 = u_in[begin]
    return isdimensionless(u0) ? (u_fac, u_in) : (u_fac/u0, ufactor!(u_in, inv(u0)))
end

function canonical_output!(u_fac::D, u_out::V) where {D<:AbstractDimensions, V<:AbstractVector{D}}
    u0 = u_out[begin]
    return isdimensionless(u0) ? (u_fac, u_out) : (u_fac*u0, ufactor!(u_out, inv(u0)))
end

function ufactor!(u::V, u_fac::D) where {D<:AbstractDimensions, V<:AbstractVector{D}}
    if ArrayInterface.can_setindex(u)
        u .= u .* u_fac
        return u
    else
        return convert(V, u .* u_fac)
    end
end
