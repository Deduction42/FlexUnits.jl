#Preamble (delete when finished)

include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")
include("linalg_types.jl")


const ArrayOfDims{D} = Union{QuantArrayDims{D}, SMatrix{N,M,D} where {N,M}} where D<:AbstractDimensions

#Dot products of dimensions
LinearAlgebra.dot(d1::AbstractDimensions, d2::AbstractDimensions) = d1*d2
dotinv(d1::AbstractDimensions, d2::AbstractDimensions) = d1/d2
function dotinv(d1::AbstractVector, d2::AbstractVector)
    length(d1) == length(d2) || throw(DimensionMismatch("Inputs had different lengths $((length(d1), length(d2)))"))
    return sum(dotinv, zip(d1, d2))
end


#DimsMap only needs to check the first row and column for equality
function Base.:(==)(d1::AbstractDimsMap, d2::AbstractDimsMap)
    equal_element(ii) = d1[begin-1+ii] == d2[begin-1+ii]
    size(d1) == size(d2) || return false
    return all(equal_element, 1:size(d1,1)) && all(equal_element, 2:size(d1, 2))
end

#For matrices, all elements must be checked for equality
Base.:(==)(d1::ArrayOfDims, d2::AbstractDimsMap) = size(d1) == size(d2) && all(==, zip(d1,d2))
Base.:(==)(d1::AbstractDimsMap, d2::ArrayOfDims) = size(d1) == size(d2) && all(==, zip(d1,d2))

function equaldims(u1::AbstractDimsMap, u2::AbstractDimsMap)
    if u1 == u2
        return u1
    else
        throw(DimensionError(u1, u2))
    end
end

Base.:+(d1::AbstractDimsMap) = d1
Base.:+(d1::AbstractDimsMap, d2::AbstractDimsMap) = equaldims(d1, d2)
Base.:-(d1::AbstractDimsMap) = d1
Base.:-(d1::AbstractDimsMap, d2::AbstractDimsMap) = equaldims(d1, d2)

#Multiplying generic factorizations with dense matrices of dimensions
Base.:*(d1::AbstractDimsMap, d2::AbstractDimsMap) = canonical!(DimsMap(u_in = uinput(d2), u_out = uoutput(d1).*dotinv(d2.uoutput, d1.uinput)))
Base.:*(d1::ArrayOfDims, d2::AbstractDimsMap) = canonical!(DimsMap(u_in = uinput(d2), u_out = d1*uoutput(d2)))
Base.:*(d1::AbstractDimsMap, d2::ArrayOfDims) = canonical!(DimsMap(u_in = inv.(d2'*inv.(uinput(d1))), u_out = uoutput(d1)))

#Multiplying specific factorizations with single dimensions
Base.:*(dm::DimsMap{<:AbstractDimensions}, d::AbstractDimensions)    = canonical!(DimsMap(u_in = dm.u_in, u_out = dm.u_out.*d))
Base.:*(d::AbstractDimensions, dm::DimsMap{<:AbstractDimensions})    = canonical!(DimsMap(u_in = dm.u_in, u_out = dm.u_out.*d))
Base.:*(dm::RepDimsMap{<:AbstractDimensions}, d::AbstractDimensions) = canonical!(RepDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))
Base.:*(d::AbstractDimensions, dm::RepDimsMap{<:AbstractDimensions}) = canonical!(RepDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))
Base.:*(dm::SymDimsMap{<:AbstractDimensions}, d::AbstractDimensions) = canonical!(SymDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))
Base.:*(d::AbstractDimensions, dm::SymDimsMap{<:AbstractDimensions}) = canonical!(SymDimsMap(u_in = dm.u_in, u_scale = dm.u_scale*d))


Base.inv(d::ArrayOfDims) = inv(DimsMap(d))

Base.:/(d1::Union{AbstractDimsMap,ArrayOfDims}, d2::Union{AbstractDimsMap,ArrayOfDims}) = d1*inv(d2)
Base.:\(d1::Union{AbstractDimsMap,ArrayOfDims}, d2::Union{AbstractDimsMap,ArrayOfDims}) = inv(d1)*d2


#======================================================================================================================
Linear algebra relationships with "Matrix" 
The outer type specializes first, so something like
    Base.inv(q::AbstractMatrix{<:Quantity}) = inv(QuantTransform(DimTransform, q))
Will be skipped in the case of Matrix{<:Quantity} (it will use inv(m::Matrix) in Base instead)
Because of this, such code will have to specify all desired concrete matrix types
======================================================================================================================#
Base.inv(q::Matrix{<:Quantity}) = inv(LinmapQuant(DimsMap, q))




#======================================================================================================================
Shortcut multiplication strategies
*(U::DimTransform, M::AbstractMatrix) = U.u_out * (U.u_in'*M)
*(M::AbstractMatrix, U::DimTransform) = (M*U.u_out) * U.u_in'
*(U1::DimTransform, U2::DimTransform) = U1.u_out.*dot(U1.u_in, U2.u_out).*U2.u_in'
      = DimTransform(u_out=U1.u_out.*dot(U1.u_in, U2.u_out), u_in=U2.u_in) 
*(U1::ScaleDimTransform, U2::ScaleDimTransform) = DimTransform(u_scale=U1.u_scale*U2.u_scale, u_out=U1.u_out) iff U1.u_out == U2.u_out
======================================================================================================================#
using Test
import .UnitRegistry.@u_str
import .UnitRegistry.@ud_str

using StaticArrays
import Random
using Statistics

#Nonlinear map
@kwdef struct PumpInput{T} <: FieldVector{2,T}
    current :: T 
    voltage :: T
end

@kwdef struct PumpOutput{T} <: FieldVector{3,T}
    power :: T 
    pressure :: T
    flow :: T 
end

function pumpfunc(x::PumpInput)
    p = x.current*x.voltage*0.9   
    return PumpOutput(power = p, pressure = sqrt(p), flow = sqrt(p))
end
pumpfunc(x::AbstractVector) = pumpfunc(PumpInput(x))


@testset "Linear Mapping Basics" begin
    Random.seed!(1234)

    #Quck tests 
    u1 = SA[u"lbf*ft", u"kW", u"rpm"]
    u2 = SA[u"kg/s", u"m^3/hr", u"kW"]

    xm = SMatrix{3,3}(randn(3,3))
    qMraw = xm.*u2./u1'
    qM = LinmapQuant(DimsMap, qMraw)

    #Matrix inversion
    x = SVector{3}(randn(3)).*u1
    y = qM*x
    @test all(x .≈ inv(qM)*y)

    #Matrix transpose
    @test all(Matrix(y') .≈ Matrix(x'*qM'))

    #Square matrices
    Σ = cov(randn(20,3)*rand(3,3))
    x = randn(3).*u2

    #Symmetric matrix
    rS = Σ.*inv.(u2).*inv.(u2)'
    qS = LinmapQuant(SymDimsMap, rS)
    @test all(qS .≈ ubase.(rS))
    @test x'*(rS)*x ≈ x'*qS*x

    #Repeatable matrix
    rR = Σ.*u2.*inv.(u2)'
    qR = LinmapQuant(RepDimsMap, rR)
    @test all(qR .≈ ubase.(rR))
    @test all((rR*rR)*x .≈ (qR*qR)*x)

    #Nonlinear mapping
    pumpunits = UnitMap(PumpInput(current=u"A", voltage=u"V"), PumpOutput(power=u"W", pressure=u"Pa", flow=u"m^3/s"))
    upumpfunc = FunctionQuant(pumpfunc, pumpunits)
    qinput = PumpInput(current=500*u"mA", voltage=6u"V")
    @test all(upumpfunc(qinput) .≈ pumpfunc(ustrip.(uinput(pumpunits), qinput)).*uoutput(pumpunits))

end
