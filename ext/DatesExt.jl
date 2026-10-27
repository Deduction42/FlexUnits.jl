module DatesExt 

import Dates: Nanosecond, Microsecond, Millisecond, Second, Minute, Hour, Day, Week, Period
import FlexUnits: seconds_unit, Quantity, Units, tobase, AbstractDimensions, AffineTransform, AbstractDimLike, StaticDims


#Convert all periods to seconds
toseconds(::Type{NanoSecond}) = 1e-9
toseconds(::Type{MicroSecond}) = 1e-6
toseconds(::Type{MilliSecond}) = 1e-3
toseconds(::Type{Second}) = true
toseconds(::Type{Minute}) = 60
toseconds(::Type{Hour}) = 3600
toseconds(::Type{Day}) = 86400
toseconds(::Type{Day}) = 604800
toseconds(p::Period) = toseconds(typeof(p))*p.value
toseconds(p::Second) = p.value

quantity(::Type{D}, p::Period) where {D <: AbstractDimLike} = toseconds(p)*seconds_unit(D)

uconvert(u::Units{D}, p::Period) where {D <: AbstractDimLike} = uconvert(u, quantity(D, p))
uconvert(u::D, p::Period) where {D <: AbstractDimLike} = uconvert(u, quantity(D, p))


end