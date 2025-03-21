

#Complex number gaps
Base.:*(q0::UnionQuantity, n::Complex) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))
Base.:*(n::Complex, q0::UnionQuantity) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))