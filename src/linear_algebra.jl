#Preamble (delete when finished)
include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")
include("linalg_types.jl")

import ArrayInterface

#Dot products of dimensions
Base.dot(d1::AbstractDimensions, d2::AbstractDimensions) = d1*d2
dotinv(d1::AbstractDimensions, d2::AbstractDimensions) = d1/d2
function dotinv(d1::AbstractVector, d2::AbstractVector)
    length(d1) == length(d2) || throw(DimensionMismatch("Inputs had different lengths $((length(d1), length(d2)))"))
    return sum(dotinv, zip(d1, d2))
end

#Canonical form has uinput[1] = dimensionless
function canonical!(u::UnitMap)
    ui1 = u.u_in[1]

    if ArrayInterface.can_setindex(u.u_in) && ArrayInterface.can_setindex(u.u_out)
        u.u_in  .= u.u_in ./ ui1
        u.u_out .= u.uout ./ ui1
        return u
    end

    return UnitMap(
        u_in = u.in./ui1, 
        u_out = u_out./ui1
    )
end

#UnitMap only needs to check the first row and column for equality
function Base.:(==)(d1::AbstractDimMap, d2::AbstractDimMap)
    equal_element(ii) = d1[begin-1+ii] == d2[begin-1+ii]
    size(d1) == size(d2) || return false
    return all(equal_element, 1:size(d1,1)) && all(equal_element, 2:size(d1, 2))
end

#For matrices, all elements must be checked for equality
Base.:(==)(d1::QuantMatrixDims, d2::AbstractDimMap) = size(d1) == size(d2) && all(==, zip(d1,d2))
Base.:(==)(d1::AbstractDimMap, d2::QuantMatrixDims) = size(d1) == size(d2) && all(==, zip(d1,d2))

function equaldims(u1::AbstractDimMap, u2::AbstractDimMap)
    if u1 == u2
        return u1
    else
        throw(DimensionError(u1, u2))
    end
end

Base.:+(d1::AbstractDimMap) = d1
Base.:+(d1::AbstractDimMap, d2::AbstractDimMap) = equaldims(d1, d2)
Base.:-(d1::AbstractDimMap) = d1
Base.:-(d1::AbstractDimMap, d2::AbstractDimMap) = equaldims(d1, d2)

Base.:*(d1::AbstractDimMap, d2::AbstractDimMap)  = canonical!(UnitMap(u_in = uinput(d2), u_out = uoutput(d1).*dotinv(d2.uoutput, d1.uinput)))
Base.:*(d1::QuantMatrixDims, d2::AbstractDimMap) = canonical!(UnitMap(u_in = uinput(d2), u_out = d1*uoutput(d2)))
Base.:*(d1::AbstractDimMap, d2::QuantMatrixDims) = canonical!(UnitMap(u_in = inv.(d2'*inv.(uinput(d1))), u_out = uoutput(d1)))

Base.inv(d::QuantMatrixDims) = inv(UnitMap(d))

Base.:/(d1::Union{AbstractDimMap,QuantMatrixDims}, d2::Union{AbstractDimMap,QuantMatrixDims}) = d1*inv(d2)
Base.:\(d1::Union{AbstractDimMap,QuantMatrixDims}, d2::Union{AbstractDimMap,QuantMatrixDims}) = inv(d1)*d2


#======================================================================================================================
Linear algebra relationships with "Matrix" 
The outer type specializes first, so something like
    Base.inv(q::AbstractMatrix{<:Quantity}) = inv(QuantTransform(DimTransform, q))
Will be skipped in the case of Matrix{<:Quantity} (it will use inv(m::Matrix) in Base instead)
Because of this, such code will have to specify all desired concrete matrix types
======================================================================================================================#
Base.inv(q::Matrix{<:Quantity}) = inv(LinmapQuant(UnitMap, q))




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
    u1 = [u"lbf*ft", u"kW", u"rpm"]
    u2 = [u"kg/s", u"m^3/hr", u"kW"]

    xm = randn(3,3)
    qM = LinmapQuant(UnitMap, xm.*u2./u1')

    #Matrix inversion
    x = randn(3).*u1
    y = qM*x
    @test all(x .≈ inv(qM)*y)

    #Matrix transpose
    @test all(Matrix(y') .≈ Matrix(x'*qM'))

    #Square matrices
    Σ = cov(randn(20,3)*rand(3,3))
    x = randn(3).*u2

    #Symmetric matrix
    rS = Σ.*inv.(u2).*inv.(u2)'
    qS = LinmapQuant(SymUnitMap, rS)
    @test all(qS .≈ ubase.(rS))
    @test x'*(rS)*x ≈ x'*qS*x

    #Repeatable matrix
    rR = Σ.*u2.*inv.(u2)'
    qR = LinmapQuant(RepUnitMap, rR)
    @test all(qR .≈ ubase.(rR))
    @test all(rR^2*x .≈ qR^2*x)

    #Nonlinear mapping
    pumpunits = UnitMap(PumpInput(current=u"A", voltage=u"V"), PumpOutput(power=u"W", pressure=u"Pa", flow=u"m^3/s"))
    upumpfunc = FunctionQuant(pumpfunc, pumpunits)
    qinput = PumpInput(current=500*u"mA", voltage=6u"V")
    @test all(upumpfunc(qinput) .≈ pumpfunc(ustrip.(uinput(pumpunits), qinput)).*uoutput(pumpunits))

end