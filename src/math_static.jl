#=
include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")
=#


#============================================================================================================================
Static dimension ops
============================================================================================================================#
Base.:(==)(d1::AbstractDimensions, d2::StaticDims) = (d1 == dimval(d2))
Base.:(==)(d1::StaticDims, d2::AbstractDimensions) = (dimval(d1) == d2)

equaldims(arg1::StaticDims, arg2::AbstractDimensions) = (dimval(arg1) == arg2) ? arg1 : throw(DimensionError((arg1,arg2)))
equaldims(arg1::AbstractDimensions, arg2::StaticDims) = (dimval(arg1) == arg2) ? arg2 : throw(DimensionError((arg1,arg2)))
equaldims(arg1::StaticDims, arg2::StaticDims) = (dimval(arg1) == dimval(arg2))  ? arg1 : throw(DimensionError((arg1,arg2)))

Base.:*(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1*D2}()
Base.:/(arg1::StaticDims{D1}, arg2::StaticDims{D2}) where {D1,D2} = StaticDims{D1/D2}()
Base.inv(arg::StaticDims{D}) where D = StaticDims{inv(D)}()
Base.:^(d::StaticDims{D}, p::Real) where D = StaticDims{D^p}()
Base.sqrt(d::StaticDims{D}) where D = StaticDims{sqrt(D)}()
Base.cbrt(d::StaticDims{D}) where D = StaticDims{sqrt(D)}()
Base.abs2(d::StaticDims{D}) where D = StaticDims{abs2(D)}()
Base.adjoint(d::StaticDims{D}) where D = StaticDims{adjoint(D)}()

@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{0}) = StaticDims{dimtype(d)()}
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{1}) = d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{2}) = d*d 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{3}) = d*d*d
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-1}) = inv(d) 
@inline Base.literal_pow(::typeof(^), d::StaticDims, ::Val{-2}) = inv(d*d)


#=
#Tests
using Test
using BenchmarkTools
import .UnitRegistry.@u_str


@testset "Static Dimensions" begin
    @test 1*StaticUnits(u"km/hr") isa Quantity{<:Any, <:StaticDims}

    #Test promotion rules
    q1 = Quantity{Float64, StaticUnits{dimension(u"km/hr"), AffineTransform}}(5, StaticUnits(u"km/hr"))
    q2 = 20.0*StaticUnits(u"m/s")
    q3 = 10*StaticUnits(u"kg/s")

    @test eltype([q1,q1]) <: Quantity{Float64, <:StaticUnits{dimension(u"m/s")}} #StaticUnits preserved
    @test eltype([q1,q2]) <: Quantity{Float64, <:StaticDims{dimension(u"m/s")}} #StaticUnits and StaticDims promote to StaticDims if dimension is the same
    @test eltype([q1,q2,q3]) <: Quantity{Float64, <:Dimensions} #Different dimensions promote to "Dimensions"

end



#=
# This is an example where static FlexUnits VASTLY outperforms Unitful due to the design choice of sticking with dimensions
import Unitful
import DynamicQuantities

v1uni  = [1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s", 1.0*Unitful.u"m/s"]
v1dyn  = [1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s", 1.0*DynamicQuantities.u"m/s"]
v1flex = [1.0*StaticUnits(u"m/s"), 1.0*StaticUnits(u"m/s"), 1.0*StaticUnits(u"m/s")]

@btime sum(x->x^2.0, v1uni)
@btime sum(x->x^2.0, $v1dyn)
@btime sum(x->x^2.0, $v1flex)
=#
=#