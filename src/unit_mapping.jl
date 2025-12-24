const UnitOrDims{D} = Union{D, AbstractUnits{D}}
const ScalarOrVec{T} = Union{T, AbstractVector{T}}

abstract type AbstractUnitMap{U<:UnitOrDims{<:AbstractDimensions}, V::ScalorOrVec{U}} <: AbstractMatrix{U} end

@kwdef struct UnitMap{U<:UnitOrDims, V<:ScalarOrVec{U}} <: AbstractUnitMap{U,V}
    u_out :: V
    u_in :: V
end

Base.getindex(m::UnitMap, ii::Integer, jj::Integer) = m.u_out[ii]/m.u_in[jj]
Base.size(m::UnitMap) = (length(m.u_out), length(m.u_in))

#Units of an adjoint are the same as an inverse, adjoint is broader because it doesn't imply one-to-one mapping
Base.inv(m::UnitMap) = UnitMap(u_out=inv.(m.u_in), u_in=inv.(m.u_out))
Base.adjoint(m::UnitMap) = inv(m)

function UnitMap(mq::AbstractMatrix{<:Union{AbstractQuantity,AbstractUnitLike}})
    u_out = dimension.(mq[:,begin])
    u_in = u_out[begin]./dimension.(mq[begin,:])

    #Check for dimensional consistency
    for jj in axes(mq,2), ii in axes(mq,1)
        (u_out[ii]/dimension(mq[ii,jj]) == u_in[jj]) || error("Unit inconsistency around index $([ii,jj]) of original matrix, expected dimension '$(u_out[ii]*inv(u_in[jj]))', found unit '$(unit(mq[ii,jj]))'")
    end

    return UnitMap(u_out=u_out, u_in=u_in)
end


@kwdef struct RUnitMap{U<:UnitOrDims, V<:ScalarOrVec{U}} <: AbstractUnitMap{U,V}
    u_scale :: U
    u_in :: V
end

Base.getindex(m::RUnitMap, ii::Integer, jj::Integer) = m.u_scale*m.u_in[ii]/m.u_in[jj]
Base.size(m::RUnitMap) = (length(m.u_in), length(m.u_in))
Base.inv(m::RUnitMap)  = RUnitMap(u_scale=inv(m.u_scale), u_in=m.u_in)

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

struct QuantTransform{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:UnitMaps{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M
    units :: U
end

function QuantTransform(::Type{U}, m::AbstractMatrix{<:Quantity}) where U <: AbstractUnitMap
    values = ustrip_base.(m)
    units  = U(dimension.(m))
    return QuantTransform(values, units)
end

Base.getindex(q::QuantTransform, ii::Integer, jj::Integer) = q.values[ii,jj] * q.units[ii,jj]
Base.size(q::QuantTransform) = size(q.values)
Base.inv(q::QuantTransform) = QuantTransform(inv(q.values), inv(q.units))


#======================================================================================================================
Linear algebra relationships with "Matrix" 
The outer type specializes first, so something like
    Base.inv(q::AbstractMatrix{<:Quantity}) = inv(QuantTransform(DimTransform, q))
Will be skipped in the case of Matrix{<:Quantity} (it will use inv(m::Matrix) in Base instead)
Because of this, such code will have to specify all desired concrete matrix types
======================================================================================================================#
Base.inv(q::Matrix{<:Quantity}) = inv(QuantTransform(UnitMap, q))


#======================================================================================================================
Shortcut multiplication strategies
*(U::DimTransform, M::AbstractMatrix) = U.u_out * (U.u_in'*M)
*(M::AbstractMatrix, U::DimTransform) = (M*U.u_out) * U.u_in'
*(U1::DimTransform, U2::DimTransform) = U1.u_out.*dot(U1.u_in, U2.u_out).*U2.u_in'
      = DimTransform(u_out=U1.u_out.*dot(U1.u_in, U2.u_out), u_in=U2.u_in) 
*(U1::ScaleDimTransform, U2::ScaleDimTransform) = DimTransform(u_scale=U1.u_scale*U2.u_scale, u_out=U1.u_out) iff U1.u_out == U2.u_out
======================================================================================================================#
