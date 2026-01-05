using StaticArrays
using BenchmarkTools

import Unitful as UF
import DynamicQuantities as DQ
using FlexUnits, .UnitRegistry

const N_ITER = Ref(10)

@kwdef struct PengRobinson{NT} <: FieldVector{10,NT}
    T  :: NT 
    P  :: NT 
    V  :: NT
    a  :: NT
    b  :: NT
    ω  :: NT
    Tc :: NT 
    Pc :: NT 
    Mw :: NT
    R  :: NT
end

PengRobinson(nt::NamedTuple) = PengRobinson(map(fn->getproperty(nt,fn), fieldnames(PengRobinson))...)

function pressure(state)
    (T, V) = (state.T, state.V)

    R = state.R
    (a, b) = (state.a, state.b)
    α  = f_alpha(state)
    r2 = sqrt(2)
    Δ  = (-1+r2, -1-r2)
    P  = R*T/(V-b) - (α*a)/((V-Δ[1]*b)*(V-Δ[2]*b))
    return P
end

function volume(state)
    (P, T) = (state.P, state.T)
    R = state.R
    V = R*T/P

    #Use the residual error of the ideal gas law to predict V and iterate
    for ii in 1:N_ITER[]
        Ph = pressure(state)
        Zr = Ph/P #Residual compressibility factor
        V  = V*Zr #Use compressibility to predict volume at P
    end

    return V 
end

function volume_function_barrier(x::PengRobinson{<:Quantity})
    x_static = (
        T=x.T|>u"K", P=x.P|>u"Pa", V=x.V|>u"m^3/mol", 
        a=x.a|>u"J*m^3/mol^2", b=x.b|>u"m^3/mol", ω=x.ω|>u"m/m", 
        Tc=x.Tc|>u"K", Pc=x.Pc|>u"Pa", Mw=x.Mw|>u"kg/mol",
        R=x.R|>u"J/(K*mol)"
    )
    return volume(x_static)
end

function f_alpha(state)
    θ  = (0.37464, 1.54226, -0.26992)
    κ  = evalpoly(state.ω, θ)
    Tr = abs(state.T)/state.Tc 
    α  = (1 + κ*(1-sqrt(Tr)))^2
    return α
end

x = PengRobinson(
    T=298, P=500e3, V=0.1, a=4.13883, b=0.000148929, ω=0.396, Tc=568.7, Pc=2.47e6, Mw=114.23e-3, R=8.31446261815324
)

#===================================================================================================
# Statically inferrable dimensions
===================================================================================================#

uf_t = (
    T=x.T*UF.u"K", P=x.P*UF.u"Pa", V=x.V*UF.u"m^3/mol", 
    a=x.a*UF.u"J*m^3/mol^2", b=x.b*UF.u"m^3/mol", ω=x.ω*UF.u"m/m", 
    Tc=x.Tc*UF.u"K", Pc=x.Pc*UF.u"Pa", Mw=x.Mw*UF.u"kg/mol",
    R=x.R*UF.u"J/(K*mol)"
)

dq_t = (
    T=x.T*DQ.u"K", P=x.P*DQ.u"Pa", V=x.V*DQ.u"m^3/mol", 
    a=x.a*DQ.u"J*m^3/mol^2", b=x.b*DQ.u"m^3/mol", ω=x.ω*DQ.u"m/m", 
    Tc=x.Tc*DQ.u"K", Pc=x.Pc*DQ.u"Pa", Mw=x.Mw*DQ.u"kg/mol",
    R=x.R*DQ.u"J/(K*mol)"
)

fl_t = (
    T=x.T*u"K", P=x.P*u"Pa", V=x.V*u"m^3/mol", 
    a=x.a*u"J*m^3/mol^2", b=x.b*u"m^3/mol", ω=x.ω*u"m/m", 
    Tc=x.Tc*u"K", Pc=x.Pc*u"Pa", Mw=x.Mw*u"kg/mol",
    R=x.R*u"J/(K*mol)"
)

println("S1.1) Iterative Peng-Robinson with statically-inferrable units\n")

print("No Units (Baseline)\t")
@btime volume($x)

print("Static Unitful.jl\t")
@btime volume($uf_t)

print("Static DynamicQ.jl\t")
@btime volume($dq_t)

print("Static FlexUnits.jl\t")
@btime volume($fl_t)

#===================================================================================================
# Dynamic dimensions
===================================================================================================#
uf_v = PengRobinson(uf_t)
dq_v = PengRobinson(dq_t)
fl_v = PengRobinson(fl_t)

println("S1.2) Iterative Peng-Robinson with dynamic units\n")

print("Dynamic Unitful  \t")
@btime volume($uf_v)

print("Dynamic DynamicQ\t")
@btime volume($dq_v)

print("Dynamic FlexUnits\t")
@btime volume($fl_v)

print("Semi-Dyn FlexUnits\t")
@btime volume_function_barrier($fl_v)