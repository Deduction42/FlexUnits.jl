using BenchmarkTools
using Random
Random.seed!(1)

# ---- Packages ----
import Unitful                            # u"..."
import DynamicQuantities                  # u"...", ua"..."
using DynamicQuantities: QuantityArray, uexpand, ustrip
using FlexUnits, .UnitRegistry            # UnitRegistry.u"..." 

println("\n== Versions ==")
@show Base.VERSION
println()

const U_R = 8.314 * Unitful.u"J/(mol*K)"
const D_R = 8.314 * DynamicQuantities.u"J/(mol*K)"
const F_R = 8.314 * UnitRegistry.u"J/(mol*K)"

# Sizes
const N  = 20000
const Ns = 1000

# ========== S1. Scalar ops (units inferable) ==========
println("S1) Scalar operations where units are inferable\n")
f1(x, y) = x * y + y * x            # result dimension is (x*y)

print("Unitful:\t")
@btime f1($(1.23 * Unitful.u"m/s"), $(0.7 * Unitful.u"m/s"));

print("DynamicQ:\t")
@btime f1($(1.23 * DynamicQuantities.u"m/s"), $(0.7 * DynamicQuantities.u"m/s"));

print("FlexU:     \t")
@btime f1($(1.23 * UnitRegistry.u"m/s"), $(0.7 * UnitRegistry.u"m/s"));

# ========== S2. Scalar ops (NON-inferable units / mixed dims) ==========
println("\nS2) Scalar ops on heterogeneous units (non-inferable)\n")

v_uni  = [1.0 * Unitful.u"m/s", 1.0 * Unitful.u"J/kg", 1.0 * Unitful.u"A/V"]
v_dyn  = [1.0 * DynamicQuantities.u"m/s", 1.0 * DynamicQuantities.u"J/kg", 1.0 * DynamicQuantities.u"A/V"]
v_flex = ubase.([1.0 * UnitRegistry.u"m/s", 1.0 * UnitRegistry.u"J/kg", 1.0 * UnitRegistry.u"A/V"])   # normalize to SI dims

print("Unitful:\t")
@btime sum(x -> x^0.0, $v_uni);
print("DynamicQ:\t")
@btime sum(x -> x^0.0, $v_dyn);
print("FlexU ubase:\t")
@btime sum(x -> x^0.0, $v_flex);

# ========== S3. Broadcasting on large arrays ==========
println("\nS3) Broadcasting on large arrays\n")

xu = randn(N) .* Unitful.u"km/s"
yu = (0.3 .+ randn(N)) .* Unitful.u"km/s"

xd = randn(N) .* DynamicQuantities.u"km/s"
yd = (0.3 .+ randn(N)) .* DynamicQuantities.u"km/s"
xqd = QuantityArray(xd);
yqd = QuantityArray(yd);
xsd = DynamicQuantities.uconvert.(DynamicQuantities.us"km/s", xd)
ysd = DynamicQuantities.uconvert.(DynamicQuantities.us"km/s", yd)

xf = randn(N) .* UnitRegistry.u"km/s"
yf = (0.3 .+ randn(N)) .* UnitRegistry.u"km/s"

xfb = ubase.(randn(N) .* UnitRegistry.u"km/s")
yfb = ubase.((0.3 .+ randn(N)) .* UnitRegistry.u"km/s")

g(x, y) = (x.^2) .+ 2.0.*x.*y .+ (y.^2)

print("Unitful:\t")
@btime g($xu, $yu);

print("DynamicQ Sym:\t")
@btime g($xsd, $ysd);

print("DynamicQ Dim:\t")
@btime g($xd, $yd);

print("DynamicQ Array:\t")
@btime g($xqd, $yqd);

print("FlexU Affine:\t")
@btime g($xf, $yf);

print("FlexU Dim:\t")
@btime g($xfb, $yfb);

# ========== S4.1. upreferred ==========
println("\nS4.1) upreferred\n")

base_array = rand(Ns)
l_uni = base_array .* Unitful.u"cm"
l_dyn = base_array .* DynamicQuantities.us"cm"
l_flex = base_array .* UnitRegistry.u"cm"

print("Unitful:\t")
@btime Unitful.upreferred.($l_uni);
print("DynamicQ:\t")
@btime uexpand.($l_dyn);
print("FlexU:  \t")
@btime ubase.($l_flex);

# ========== S4.2. ustrip ==========
println("\nS4.2) ustrip\n")
print("Unitful:\t")
@btime Unitful.ustrip.(Unitful.u"mm", $l_uni);
print("DynamicQ:\t")
@btime DynamicQuantities.ustrip.(DynamicQuantities.u"mm", $l_dyn);
print("FlexU:  \t")
@btime FlexUnits.ustrip.(UnitRegistry.u"mm", $l_flex);

# ========== S4.3. uconvert ==========
println("\nS4.3) uconvert to arbitrary units\n")
print("Unitful:\t")
@btime Unitful.uconvert.(Unitful.u"ft", $l_uni);
print("DynamicQ:\t")
@btime DynamicQuantities.uconvert.(DynamicQuantities.us"ft", $l_dyn);
print("FlexU:  \t")
@btime FlexUnits.uconvert.(UnitRegistry.u"ft", $l_flex);

# ========== S5. Affine units (°C/°F) ==========
println("\nS5) Affine units (°C/°F) handling: PV = nRT at 25°C, 101.3kPa, n=1 mol\n")
Tc = 25 .+ randn(Ns);
Tu = (Tc .+ 273.15) .* Unitful.u"K";
Td = [Tci * DynamicQuantities.ua"degC" for Tci in Tc];
Tf = ubase.(Tc .* UnitRegistry.u"°C");

pu, nu = 101.3 * Unitful.u"kPa", 1.0 * Unitful.u"mol";
pd, nd = 101.3 * DynamicQuantities.u"kPa", 1.0 * DynamicQuantities.u"mol";
pf, nf = ubase(101.3 * UnitRegistry.u"kPa"), ubase(1.0 * UnitRegistry.u"mol");

print("Unitful:\t")
@btime ($nu .* U_R .* $Tu) ./ $pu;

print("DynamicQ:\t")
@btime ($nd .* D_R .* $Td) ./ $pd;

print("FlexU  :\t")
@btime ($nf .* F_R .* $Tf) ./ $pf;

# ========== S6. Struct storage & field access ==========
println("\nS6) Structs with quantities as fields")

# Unitful: keep fields concrete via type params (typical pattern)
const UQK   = typeof(1.0 * Unitful.u"K")
const UQkPa = typeof(1.0 * Unitful.u"kPa")
struct UState{TK,TP}
    T::TK
    p::TP
end

# DynamicQuantities single concrete type works
const UD = typeof(1.0*DynamicQuantities.u"K")
struct DState
    T::UD
    p::UD
end

# FlexUnits single concrete type works
const UF = typeof(ubase(1.0*UnitRegistry.u"K"))
struct FState
    T::UF
    p::UF
end

function build_states_unitful(n)
    v = Vector{UState{UQK,UQkPa}}(undef, n)
    @inbounds for i = 1:n
        v[i] = UState((290.0 + rand() * 20) * Unitful.u"K", (90.0 + rand() * 30) * Unitful.u"kPa")
    end
    v
end

function build_states_dyn(n)
    [DState((290.0 + rand() * 20) * DynamicQuantities.u"K",
        (90.0 + rand() * 30) * DynamicQuantities.u"kPa") for _ = 1:n]
end

function build_states_flex(n)
    [FState((290.0 + rand() * 20) * UnitRegistry.u"K",
        (90.0 + rand() * 30) * UnitRegistry.u"kPa") for _ = 1:n]
end

println("\nS6.1) Construct\n")
print("Unitful:\t")
@btime build_states_unitful($Ns);
print("DynamicQ:\t")
@btime build_states_dyn($Ns);
print("FlexU:  \t")
@btime build_states_flex($Ns);

println("\nS6.2) Access\n")
sumT_u(v) = sum(s -> s.T, v)
sumT_d(v) = sum(s -> s.T, v)
sumT_f(v) = sum(s -> s.T, v)

vu = build_states_unitful(Ns);
vd = build_states_dyn(Ns);
vf = build_states_flex(Ns);

print("Unitful:\t")
@btime sumT_u($vu);
print("DynamicQ:\t")
@btime sumT_d($vd);
print("FlexU:  \t")
@btime sumT_f($vf);

println("\nS6.3) Heterogeneous states (different units across elements)\n")
vdh = [DState(300 * DynamicQuantities.u"K", 100 * DynamicQuantities.u"kPa"),
    DState(27 * DynamicQuantities.ua"degC", 14 * DynamicQuantities.u"Pa")]
vfh = [FState(300 * UnitRegistry.u"K", 100 * UnitRegistry.u"kPa"),
    FState(27 * UnitRegistry.u"°C", 14 * UnitRegistry.u"Pa")]
vuh = Any[
    UState(300 * Unitful.u"K", 100 * Unitful.u"kPa"),
    UState((273.15 + 27) * Unitful.u"K", 14 * Unitful.u"Pa")
]
print("Unitful:\t")
@btime sumT_u($vuh);
print("DynamicQ:\t")
@btime sumT_d($vdh);
print("FlexU:  \t")
@btime sumT_f($vfh);