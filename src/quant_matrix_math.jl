#Preamble (delete when finished)
include("fixed_rational.jl")
include("types.jl")
include("utils.jl")
include("conversions.jl")
include("math.jl")
include("RegistryTools.jl")
include("UnitRegistry.jl")
include("quant_matrix_types.jl")





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