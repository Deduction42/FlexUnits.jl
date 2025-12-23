const UnitOrDims{D} = Union{D, AbstractUnits{D}}

abstract type AbstractUnitMap{U<:UnitOrDims{<:AbstractDimensions}} <: AbstractMatrix{U} end

@kwdef struct UnitMap{U<:AbstractUnitLike, V<:AbstractVector{U}} <: AbstractUnitMap{U}
    u_out :: V
    u_in :: V
end

Base.getindex(m::UnitMap, ii::Integer, jj::Integer) = m.u_out[ii]*m.u_in[jj]
Base.size(m::UnitMap) = (length(m.u_out), length(m.u_in))

#Units of an adjoint are the same as an inverse, adjoint is broader because it doesn't imply one-to-one mapping
Base.adjoint(m::UnitMap) = UnitMap(u_out=inv.(m.u_in), u_in=inv.(m.u_out))
Base.inv(m::UnitMap) = adjoint(m)

function UnitMap(mq::AbstractMatrix{<:Union{AbstractQuantity,AbstractUnitLike}})
    u_out = dimension.(mq[:,begin])
    u_in = dimension.(mq[begin,:])./u_out[begin]

    #Check for dimensional consistency
    for jj in axes(mq,2), ii in axes(mq,1)
        (dimension(mq[ii,jj])/u_out[ii] == u_in[jj]) || error("Unit inconsistency around index $([ii,jj]), expected dimension of $(u_out[ii]*u_in[jj])")
    end

    return UnitMap(u_out=u_out, u_in=u_in)
end


@kwdef struct RUnitMap{D<:AbstractDimensions, V<:AbstractVector{D}} <: AbstractUnitMap{D}
    uscale :: D
    u_in :: V
end

Base.getindex(m::RUnitMap, ii::Integer, jj::Integer) = m.uscale*m.u_in[ii]/m.u_in[jj]
Base.size(m::RUnitMap) = (length(m.u_in), length(m.u_in))
Base.inv(m::RUnitMap)  = RUnitMap(uscale=inv(m.uscale), u_in=m.u_in)

function RUnitMap(md::UnitMap)
    #Matrix must be square
    sz = size(md)
    sz[1] == sz[2] || throw(DimensionMismatch("Recursive Unit Mapping must be square: dimensions are $(sz)"))

    #Calculate the uniform scale
    uscale = md.u_out[begin]*md.u_in[begin]

    #Verify that uniform scale is consistent 
    for (u_out, u_in) in zip(md.u_out, md.u_in)
        u_out*u_in == uscale || error("Cannot convert to Recursive Unit Mapping: $(md.u_out) and $(md.u_in) must share a common factor")
    end

    return RUnitMap(uscale=uscale, u_in=md.u_out./uscale)
end

function RUnitMap(mq::AbstractMatrix{<:Union{AbstractQuantity,AbstractUnitLike}})
    return RUnitMap(UnitMap(mq))
end

UnitMap(md::RUnitMap) = UnitMap(u_out=md.u_in.*md.uscale, u_in=inv.(md.u_in))

struct QuantTransform{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:AbstractUnitMap{D}} <: AbstractMatrix{Quantity{T,D}}
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
*(U1::ScaleDimTransform, U2::ScaleDimTransform) = DimTransform(uscale=U1.uscale*U2.uscale, u_out=U1.u_out) iff U1.u_out == U2.u_out
======================================================================================================================#
