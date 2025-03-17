"""
    AffineTransformation(target::AbstractUnitLike, current::AbstractUnitLike)

Calculates the affine transformation to convert a value from the current units to the target units
"""
@kwdef struct AffineTransformation
    scale :: Float64
    offset :: Float64
end

function AffineTransformation(target::AbstractUnitLike, current::AbstractUnitLike)
    return AffineTransformation(
        scale  = uscale(current)/uscale(target),
        offset = (uoffset(current) - uoffset(target))/uscale(target)
    )
end

Base.muladd(x, affine::AffineTransformation) = muladd(x, affine.scale, affine.offset)

convert(::Type{<:Quantity}, q::UnionQuantity) = Quantity(q)
convert(::Type{<:NumberQuantity}, q::UnionQuantity) = NumberQuantity(q)
convert(::Type{<:RealQuantity}, q::UnionQuantity) = RealQuantity(q)



#=
uconvert(u::Units, q::UnionQuantity) will calculate the affine transformation, apply it and return the narrowest quantity in new units

convert(::Quantity{Dimensions}, q::UnionQuantity) will uconvert(dimension(q), q) then apply Quantity
convert(::Quantity{ScalarUnits}, q::UnionQuantity) will unconvert(ScalarUnits(scale=scale(q), dims=dimmension(q)), q), then apply Quantity
convert(::Quantity{AffineUnits}, q::UnionQuantity) will ....
convert(::Quantity{Dimensions}, q::AffineUnits) will fail if offset is not zero

can potentially define:
convert(::Type{Q}, q::UnionQuantity) where {U<:Dimensions,Q<:UnionQuantity{U}} = Q(uconvert(dimension(U), q)
=#
