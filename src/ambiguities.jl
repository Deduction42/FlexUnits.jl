#Complex number gaps for *
Base.:*(q0::UnionQuantity, n::Complex) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))
Base.:*(n::Complex, q0::UnionQuantity) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))
Base.:*(q0::UnionQuantity, n::Complex{Bool}) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))
Base.:*(n::Complex{Bool}, q0::UnionQuantity) = (q = ubase(q0); quantity(ustrip(q)*n, unit(q)))

#Complex number gaps for /
Base.:/(q0::UnionQuantity, n::Complex) = (q = ubase(q0); quantity(ustrip(q)/n, unit(q)))
Base.:/(n::Complex, q0::UnionQuantity) = (q = ubase(q0); quantity(n/ustrip(q), inv(unit(q))))
Base.:/(q0::UnionQuantity, n::Complex{Bool}) = (q = ubase(q0); quantity(ustrip(q)/n, unit(q)))
Base.:/(n::Complex{Bool}, q0::UnionQuantity) = (q = ubase(q0); quantity(n/ustrip(q), inv(unit(q))))

#Complex number gaps for +
Base.:+(q0::UnionQuantity, n::Complex) = dimensionless(q0) + n 
Base.:+(n::Complex, q0::UnionQuantity) = dimensionless(q0) + n 
Base.:+(q0::UnionQuantity, n::Complex{Bool}) = dimensionless(q0) + n 
Base.:+(n::Complex{Bool}, q0::UnionQuantity) = dimensionless(q0) + n

#Complex number gaps for -
Base.:-(q0::UnionQuantity, n::Complex) = dimensionless(q0) - n 
Base.:-(n::Complex, q0::UnionQuantity) = n - dimensionless(q0) 
Base.:-(q0::UnionQuantity, n::Complex{Bool}) = dimensionless(q0) - n 
Base.:-(n::Complex{Bool}, q0::UnionQuantity) = n - dimensionless(q0)

