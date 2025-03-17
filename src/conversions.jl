@kwdef struct AffineTransformation
    scale :: Float64
    offset :: Float64
end

Base.muladd(x, affine::AffineTransformation) = muladd(x, affine.scale, affine.offset)

"""
    AffineTransformation(target::AbstractUnitLike, current::AbstractUnitLike)

Calculates the affine transformation to convert a value from the current units to the target units
"""
function AffineTransformation(target::AbstractUnitLike, current::AbstractUnitLike)
    return AffineTransformation(
        scale  = uscale(current)/uscale(target),
        offset = (offset(current) - offset(target))/scale(target)
    )
end


convert(::Type{<:Quantity}, q::UnionQuantity) = Quantity(q)
convert(::Type{<:NumberQuantity}, q::UnionQuantity) = NumberQuantity(q)
convert(::Type{<:RealQuantity}, q::UnionQuantity) = RealQuantity(q)

