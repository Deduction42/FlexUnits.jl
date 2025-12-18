abstract type AbstractDimTransform{D<:AbstractDimensions} <: AbstractMatrix{D} end

@kwdef struct DimTransform{D<:AbstractDimensions, V<:AbstractVector{D}} <: AbstractDimTransform{D}
    colfac :: V
    rowfac :: V
end

Base.getindex(m::DimTransform, ii::Integer, jj::Integer) = m.colfac[ii]*m.rowfac[jj]
Base.size(m::DimTransform) = (length(m.colfac), length(m.rowfac))

#Units of an adjoint are the same as an inverse, adjoint is broader because it doesn't imply one-to-one mapping
Base.adjoint(m::DimTransform) = DimTransform(colfac=inv.(m.rowfac), rowfac=inv.(m.colfac))
Base.inv(m::DimTransform) = adjoint(m)

function DimTransform(mq::AbstractMatrix{<:Union{AbstractQuantity,AbstractUnitLike}})
    colfac = dimension.(mq[:,begin])
    rowfac = dimension.(mq[begin,:])./colfac[begin]

    #Check for dimensional consistency
    for jj in axes(mq,2), ii in axes(mq,1)
        (dimension(mq[ii,jj])/colfac[ii] == rowfac[jj]) || error("Unit inconsistency around index $([ii,jj]), expected dimension of $(colfac[ii]*rowfac[jj])")
    end

    return DimTransform(colfac=colfac, rowfac=rowfac)
end


@kwdef struct ScaleDimTransform{D<:AbstractDimensions, V<:AbstractVector{D}} <: AbstractDimTransform{D}
    uscale :: D
    colfac :: V
end

Base.getindex(m::ScaleDimTransform, ii::Integer, jj::Integer) = m.uscale*m.colfac[ii]/m.colfac[jj]
Base.size(m::ScaleDimTransform) = (length(m.colfac), length(m.colfac))
Base.inv(m::ScaleDimTransform)  = ScaleDimTransform(uscale=inv(m.uscale), colfac=m.colfac)

function ScaleDimTransform(md::DimTransform)
    #Matrix must be square
    sz = size(md)
    sz[1] == sz[2] || throw(DimensionMismatch("matrix is not square: dimensions are $(sz)"))

    #Calculate the uniform scale
    uscale = md.colfac[begin]*md.rowfac[begin]

    #Verify that uniform scale is consistent 
    for (colfac, rowfac) in zip(md.colfac, md.rowfac)
        colfac*rowfac == uscale || error("Cannot convert DimTransform to ScaleDimTransform, $(md.colfac) and $(md.rowfac) don't multiply to a common unit")
    end

    return ScaleDimTransform(uscale=uscale, colfac=md.colfac./uscale)
end

function ScaleDimTransform(mq::AbstractMatrix{<:Union{AbstractQuantity,AbstractUnitLike}})
    return ScaleDimTransform(DimTransform(mq))
end

DimTransform(md::ScaleDimTransform) = DimTransform(colfac=md.colfac.*md.uscale, rowfac=inv.(md.colfac))

struct QuantTransform{T, D<:AbstractDimensions, M<:AbstractMatrix{T}, U<:AbstractDimTransform{D}} <: AbstractMatrix{Quantity{T,D}}
    values :: M 
    units :: U
end

function QuantTransform(::Type{U}, m::AbstractMatrix{<:Quantity}) where U <: AbstractDimTransform
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
Base.inv(q::Matrix{<:Quantity}) = inv(QuantTransform(DimTransform, q))


#======================================================================================================================
Shortcut multiplication strategies
*(U::DimTransform, M::AbstractMatrix) = U.colfac * (U.rowfac'*M)
*(M::AbstractMatrix, U::DimTransform) = (M*U.colfac) * U.rowfac'
*(U1::DimTransform, U2::DimTransform) = U1.colfac.*dot(U1.rowfac, U2.colfac).*U2.rowfac'
      = DimTransform(colfac=U1.colfac.*dot(U1.rowfac, U2.colfac), rowfac=U2.rowfac) 
*(U1::ScaleDimTransform, U2::ScaleDimTransform) = DimTransform(uscale=U1.uscale*U2.uscale, colfac=U1.colfac) iff U1.colfac == U2.colfac
======================================================================================================================#
