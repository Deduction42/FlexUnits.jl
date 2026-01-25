module ForwardDiffExt
import FlexUnits
import ForwardDiff

function Base.convert(d::Type{ForwardDiff.Dual{T, V, N}}, q::FlexUnits.Quantity{<:Real, <:FlexUnits.AbstractDimLike}) where {T, V, N}
    return d(FlexUnits.dstrip(q))
end

ForwardDiff.can_dual(::Type{Q}) where {Q<:FlexUnits.Quantity{<:Real, <:FlexUnits.AbstractDimLike}} = true

end