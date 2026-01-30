using LinearAlgebra
using StaticArrays
using SparseArrays

import ArrayInterface
import StaticArrays.StaticLUMatrix
import SparseArrays.AbstractSparseMatrixCSC
import SparseArrays.AbstractCompressedVector

# AbstractDimsMap is similar to AbstractArray, but doesn't subtype to it; it subtypes to AbstractUnitMap which is not really an array
# Moreover subtyping to AbstractArray produces ambiguity problems, it's easier to just redefine the API methods
abstract type AbstractDimsMap{D<:AbstractDimensions} <: AbstractUnitMap{D} end

Base.eltype(::Type{<:AbstractDimsMap{D}}) where D = D
Base.IndexStyle(::Type{AbstractDimsMap}) = IndexCartesian()
Base.getindex(m::AbstractDimsMap, ind::CartesianIndex{2}) = m[ind[1],ind[2]]
Base.getindex(m::AbstractDimsMap, ind::Integer) = m[CartesianIndices(m)[ind]]
Base.CartesianIndices(m::AbstractDimsMap) = CartesianIndices(axes(m))
Base.length(m::AbstractDimsMap) = prod(size(m))
Base.collect(m::AbstractDimsMap) = uoutput(m).*inv.(uinput(m)')
Base.axes(m::AbstractDimsMap, d::Integer) = d <= 2 ? axes(m)[d] : OneTo(1)




"""
    AdjointDmap{D, M<:AbstractUnitMap{D}} <: AbstractUnitMap{D}

Wraps a unitmap as an adjoint/transpose (and avoids subtyping to AbstractArray)
"""
struct AdjointDmap{D, M<:AbstractDimsMap{D}} <: AbstractDimsMap{D}
    parent :: M 
end
Base.IndexStyle(::Type{AdjointDmap{D,M}}) where {D,M} = IndexStype(M)
Base.transpose(m::AbstractUnitMap) = AdjointDmap(m)
Base.transpose(m::AdjointDmap) = m.parent
Base.adjoint(m::AbstractUnitMap) = AdjointDmap(m)
Base.adjoint(m::AdjointDmap) = m.parent 
Base.getindex(m::AdjointDmap, ind1::Integer, ind2::Integer) = getindex(m.parent, ind2, ind1)
Base.axes(m::AdjointDmap) = reverse(axes(m.parent))
uoutput(m::AdjointDmap) = inv.(uinput(m.parent))*uscale(m.parent)
uoutput(m::AdjointDmap) = inv.(uinput(m.parent))*uscale(m.parent)
uinput(m::AdjointDmap) = inv.(uoutput(m.parent))*uscale(m.parent)
uscale(m::AdjointDmap) = uscale(m.parent)

#Base.iterate(m::AbstractDimsMap, i=1) = (@inline; (i - 1)%UInt < length(m)%UInt ? (m[i], i + 1) : nothing)

"""
    ArrayDims{D<:AbstractDimensions, A<:AbstractArray} <: AbstractArray{D}

Wraps a quantity matrix A so that getting its index returns the dimensiosn of the element
"""
struct QuantArrayDims{D<:AbstractDimensions, N, A<:AbstractArray} <: AbstractArray{D,N}
    array :: A
    function QuantArrayDims(a::AbstractArray{<:Any,N}) where N
        D = dimvaltype(eltype(a))
        return new{D,N,typeof(a)}(a)
    end
end
Base.IndexStyle(::Type{QuantArrayDims{D,A}}) where{D,A} = IndexStyle(A)
Base.getindex(m::QuantArrayDims, args...) = broadcast(dimension, getindex(m.array, args...))
Base.size(m::QuantArrayDims) = size(m.array)
Base.axes(m::QuantArrayDims) = axes(m.array)
ArrayInterface.can_setindex(::Type{QuantArrayDims}) = false
dimension(m::AbstractArray) = QuantArrayDims(m)
dimension(m::SArray) = dimension.(m)
dstrip(m::AbstractArray) = dstrip.(m)

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
    u_in  :: TI
    u_out :: TO
end

Used to represent a unit transformation from input dimensions 'u_in' to outpout dimensions 'u_out'.
This is like a unit map but focuses on dimensions, simplifying linear algebra (it subtypes to Matrix)
"""
@kwdef struct DimsMap{D<:AbstractDimensions, TI<:AbstractVector{D}, TO<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_scale :: D
    u_in  :: TI
    u_out :: TO
    function DimsMap{D, TI, TO}(u_scale, u_in_raw, u_out_raw) where {D<:AbstractDimensions, TI<:AbstractVector{D}, TO<:AbstractVector{D}}
        (u_scale, u_in)  = canonical_input(u_scale, u_in_raw)
        (u_scale, u_out) = canonical_output(u_scale, u_out_raw)
        return new{D, TI, TO}(u_scale, u_in, u_out)
    end
    function DimsMap(u_scale::D, u_in::TI, u_out::TO) where {D<:AbstractDimensions, TI<:AbstractVector{D}, TO<:AbstractVector{D}}
        return new{D, TI, TO}(u_scale, u_in, u_out)
    end
end

Base.axes(m::DimsMap) = (axes(m.u_out)[1], axes(m.u_in)[1])
Base.getindex(m::DimsMap, ii::Integer, jj::Integer) = m.u_out[ii]/m.u_in[jj]*m.u_scale
Base.size(m::DimsMap) = (length(m.u_out), length(m.u_in))
Base.inv(m::DimsMap) = DimsMap(u_scale=m.u_scale, u_out=m.u_in, u_in=m.u_out)
#Base.adjoint(m::DimsMap) = DimsMap(u_out=inv.(m.u_in).*m.u_out[begin], u_in=inv.(m.u_out).*m.u_out[begin])
#Base.transpose(m::DimsMap) = adjoint(m)
uoutput(m::DimsMap) = m.u_out
uinput(m::DimsMap) = m.u_in
uscale(m::DimsMap) = m.u_scale

function Base.firstindex(m::DimsMap, d) 
    if d==1 
        return firstindex(m.u_out) 
    elseif d==2 
        return firstindex(m.u_in)
    end
    return 1
end

function DimsMap(md::AbstractMatrix{<:AbstractDimensions})
    u_scale = md[begin,begin]
    u_out = md[:,begin]./u_scale
    u_in  = u_scale./md[begin,:]

    #Check for dimensional consistency
    for jj in axes(md,2), ii in axes(md,1)
        md_ij = u_out[ii]/u_in[jj]*u_scale
        (md_ij == md[ii, jj]) || error("Unit inconsistency around index $([ii, jj]) of original matrix, expected dimension '$(md_ij))', found dimension '$(md[ii, jj])'")
    end
    return DimsMap(u_scale=u_scale, u_out=u_out, u_in=u_in)
end

DimsMap(mq::AbstractMatrix{<:QuantUnion}) = DimsMap(QuantArrayDims(mq))

#=
function canonical!(u::DimsMap)
    u0 = u.u_in[begin]
    isdimensionless(u0) && return u
    
    unew = if ArrayInterface.can_setindex(u.u_in) && ArrayInterface.can_setindex(u.u_out)
        u.u_in  .= u.u_in ./ u0
        u.u_out .= u.u_out ./ u0
        u
    else
        DimsMap(
            u_in = u.u_in./u0, 
            u_out = u.u_out./u0
        )
    end
    return unew
end
=#

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
    function RepDimsMap{D, TI}(u_scale, u_in_raw) where {D<:AbstractDimensions, TI<:AbstractVector{D}}
        (u_scale, u_in)  = canonical_rep(u_scale, u_in_raw)
        return new{D, TI}(u_scale, u_in)
    end
    function DimsMap(u_scale::D, u_in::TI) where {D<:AbstractDimensions, TI<:AbstractVector{D}}
        return new{D, TI}(u_scale, u_in)
    end
end

Base.axes(m::RepDimsMap) = (axes(m.u_in)[1], axes(m.u_in)[1])
Base.getindex(m::RepDimsMap, ii::Integer, jj::Integer) = m.u_scale*m.u_in[ii]/m.u_in[jj]
Base.size(m::RepDimsMap) = (length(m.u_in), length(m.u_in))
Base.inv(m::RepDimsMap)  = RepDimsMap(u_scale=inv(m.u_scale), u_in=m.u_in)
#Base.adjoint(m::RepDimsMap) = RepDimsMap(u_scale=m.u_scale, u_in=inv.(m.u_in))
#Base.transpose(m::RepDimsMap) = adjoint(m)
Base.firstindex(m::RepDimsMap, d) = firstindex(m.u_in)
uoutput(m::RepDimsMap) = m.u_in
uinput(m::RepDimsMap) = m.u_in
uscale(m::RepDimsMap) = m.u_scale

function RepDimsMap(md::DimsMap)
    #Matrix must be square
    sz = size(md)
    sz[1] == sz[2] || throw(DimensionMismatch("Repeatable Unit Mapping must be square: dimensions are $(sz)"))

    #Ensure inputs and outputs are equal
    md.u_in == md.u_out || error("Cannot convert to Repeatable Unit Mapping: $(md.u_out) and $(md.u_in) must be equal")

    return RepDimsMap(u_scale=md.u_scale, u_in=md.u_in)
end

RepDimsMap(mq::AbstractMatrix{<:Union{QuantUnion,AbstractUnitLike}}) = RepDimsMap(DimsMap(mq))
RepDimsMap(md::RepDimsMap) = md
DimsMap(md::RepDimsMap) = DimsMap(u_out=md.u_in.*md.u_scale, u_in=md.u_in)

#=
function canonical!(u::RepDimsMap)
    u0 = u.u_in[begin]
    isdimensionless(u0) && return u

    u_in = if ArrayInterface.can_setindex(u.u_in)
        u.u_in .= u.u_in ./ u0
        u.u_in
    else
        u.in./u0
    end
    return RepDimsMap(u_in = u_in, u_scale = u.u_scale)
end
=#

"""
struct SymUnitMap{D<:AbstractDimensions, TI<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_scale :: D
    u_in :: TI
end

Used to represent a special kind of dimensional transformation (a symmetric transformation) of dimensions 'u_in'.
If "U" is symmetric, then "x'U*x" is a valid operation and the dimensions of "U*x" are similar to the inverse of the
dimensions of "x". This structure enables certain kinds of operations reserved for symmetric matrices.
"""
@kwdef struct SymDimsMap{D<:AbstractDimensions, TI<:AbstractVector{D}} <: AbstractDimsMap{D}
    u_scale :: D
    u_in :: TI
    function SymDimsMap{D, TI}(u_scale, u_in_raw) where {D<:AbstractDimensions, TI<:AbstractVector{D}}
        (u_scale, u_in)  = canonical_sym(u_scale, u_in_raw)
        return new{D, TI}(u_scale, u_in)
    end
    function SymDimsMap(u_scale::D, u_in::TI) where {D<:AbstractDimensions, TI<:AbstractVector{D}}
        return new{D, TI}(u_scale, u_in)
    end
end

Base.axes(m::SymDimsMap) = (axes(m.u_in)[1], axes(m.u_in)[1])
Base.getindex(m::SymDimsMap, ii::Integer, jj::Integer) = m.u_scale/(m.u_in[ii]*m.u_in[jj])
Base.size(m::SymDimsMap) = (length(m.u_in), length(m.u_in))
Base.inv(m::SymDimsMap)  = SymDimsMap(u_scale=inv(m.u_scale), u_in=inv.(m.u_in))
Base.adjoint(m::SymDimsMap) = m
Base.transpose(m::SymDimsMap) = m
Base.firstindex(m::SymDimsMap, d) = firstindex(m.u_in)
uoutput(m::SymDimsMap) = map(u->inv(u), m.u_in)
uinput(m::SymDimsMap) = m.u_in
uscale(m::SymDimsMap) = m.u_scale

function SymDimsMap(md::DimsMap)
    #Matrix must be square
    sz = size(md)
    sz[1] == sz[2] || throw(DimensionMismatch("Symmetric Unit Mapping must be square: dimensions are $(sz)"))

    #Check if the inverses are eaual
    md.u_out[begin] == inv.(md.u_in[begin]) || error("Cannot convert to Symmetric Unit Mapping: $(md.u_in) and $(md.u_out) must be inverses of each other")

    return SymDimsMap(u_scale=md.u_scale, u_in=md.u_in)
end

SymDimsMap(mq::AbstractMatrix{<:Union{QuantUnion,AbstractUnitLike}}) = SymDimsMap(DimsMap(mq))
SymDimsMap(md::SymDimsMap) = md
DimsMap(md::SymDimsMap) = DimsMap(u_out=md.u_in.*md.u_scale, u_in=md.u_in)

#=
function canonical!(u::SymDimsMap)
    u0 = u.u_in[begin]
    isdimensionless(u0) && return u

    u_in = if ArrayInterface.can_setindex(u.u_in)
        u.u_in .= u.u_in ./ u0
        u.u_in
    else
        u.in./u0
    end
    return SymDimsMap(u_in = u_in, u_scale = u.u_scale/u0^2)
end
=#
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
    new_m = todims.(m, u.u_out./u.u_in')
    new_u = DimsMap(u_in = dimension.(u.u_in), u_out = dimension.(u.u_out))
    return LinmapQuant(new_m, canonical!(new_u))
end

function LinmapQuant(::Type{U}, m::AbstractMatrix) where U <: AbstractDimsMap
    values = dstrip.(m)
    units  = U(dimension(m))
    return LinmapQuant(values, units)
end

LinmapQuant(m::AbstractMatrix) = LinmapQuant(DimsMap, m)
LinmapQuant(m::LinmapQuant) = m

ustrip(lq::LinmapQuant) = lq.values
dstrip(lq::LinmapQuant) = lq.values
unit(lq::LinmapQuant) = lq.dims
dimension(lq::LinmapQuant) = lq.dims
ubase(lq::LinmapQuant) = lq

Base.IndexStyle(::Type{<:LinmapQuant}) = IndexCartesian()
Base.getindex(q::LinmapQuant, ii::Integer, jj::Integer) = q.values[ii,jj] * q.dims[ii,jj]
Base.size(q::LinmapQuant) = size(q.values)
Base.inv(q::LinmapQuant) = LinmapQuant(inv(q.values), inv(q.dims))
Base.transpose(q::LinmapQuant) = LinmapQuant(transpose(q.values), transpose(q.dims))
Base.adjoint(q::LinmapQuant) = LinmapQuant(adjoint(q.values), adjoint(q.dims))



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

VectorQuant(v::AbstractVector) = VectorQuant(dstrip.(v), dimension.(v))
VectorQuant(v::VectorQuant) = v

ustrip(lq::VectorQuant) = lq.values
dstrip(lq::VectorQuant) = lq.values
unit(lq::VectorQuant) = lq.dims
dimension(lq::VectorQuant) = lq.dims
ubase(lq::VectorQuant) = lq

Base.IndexStyle(::Type{<:VectorQuant{<:Any, <:Any, V}}) where {V} = IndexStyle(V)
Base.getindex(q::VectorQuant, ii::Integer) = q.values[ii] * q.dims[ii]
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
unit(fq::FactorQuant) = getfield(fq, :dims)
dimension(fq::FactorQuant) = getfield(fq, :dims)
Base.inv(fq::FactorQuant) = LinmapQuant(inv(ustrip(fq)), inv(unit(fq)))
LinearAlgebra.inv!(fq::FactorQuant) = LinmapQuant(inv!(ustrip(fq)), inv(unit(fq)))

# LU Factorization ===================================================================================
LinearAlgebra.lu(mq::LinmapQuant; kwargs...) = FactorQuant(lu(ustrip(mq); kwargs...), unit(mq))
LinearAlgebra.lu(mq::LinmapQuant, ::Val{true}; kwargs...) = FactorQuant(lu(ustrip(mq), Val(true); kwargs...), unit(mq))
LinearAlgebra.lu(mq::LinmapQuant, ::Val{false}; kwargs...) = FactorQuant(lu(ustrip(mq), Val(false); kwargs...), unit(mq))


#May need to iterate over more subtypes of AbstractMatrix
StaticArrays.lu(mq::StaticLUMatrix{N,M,<:Quantity}; kwargs...) where {N,M} = lu(LinmapQuant(DimsMap, mq); kwargs...)
LinearAlgebra.lu(mq::AbstractMatrix{<:Quantity}; kwargs...) = lu(LinmapQuant(DimsMap, mq); kwargs...)
LinearAlgebra.lu(mq::AbstractMatrix{<:Quantity}, ::Val{true}; kwargs...) = lu(LinmapQuant(DimsMap, mq), Val(true); kwargs...)
LinearAlgebra.lu(mq::AbstractMatrix{<:Quantity}, ::Val{false}; kwargs...) = lu(LinmapQuant(DimsMap, mq), Val(false); kwargs...)


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


#======================================================================================================================
Special cases
======================================================================================================================#
#UniformScaling with dynamic dimensions should produce unknown dimension on off-diagonals (consistent with other behaviour)
Base.getindex(J::UniformScaling{T}, i::Integer, j::Integer) where T<:Quantity = ifelse(i==j, J.λ, zero(T))


#======================================================================================================================
Utility functions
======================================================================================================================#
function canonical_input(u_scale::D, u_in::V) where {D<:AbstractDimensions, V<:AbstractVector{D}}
    u0 = u_in[begin]
    return isdimensionless(u0) ? (u_scale, u_in) : (u_scale/u0, convert(V, u_in./u0))
end

function canonical_output(u_scale::D, u_out::V) where {D<:AbstractDimensions, V<:AbstractVector{D}}
    u0 = u_out[begin]
    return isdimensionless(u0) ? (u_scale, u_out) : (u_scale*u_0, convert(V, u_out./u0))
end

function canonical_rep(u_scale::D, u_in::V) where {D<:AbstractDimensions, V<:AbstractVector{D}}
    u0 = u_in[begin]
    return isdimensionless(u0) ? (u_scale, u_in) : (u_scale, convert(V, u_in./u0))
end

function canonical_sym(u_scale::D, u_in::V) where {D<:AbstractDimensions, V<:AbstractVector{D}}
    u0 = u_in[begin]
    return isdimensionless(u0) ? (u_scale, u_in) : (u_scale/u0^2, convert(V, u_in./u0))
end