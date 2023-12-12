def coron_trivariate(pol, X, Y, Z, l=2):
    P.<x,y,z> = PolynomialRing(ZZ)

    variables = [str(pol.variables()[0]), str(pol.variables()[1]), str(pol.variables()[2])]

    if ("x1" in variables) and ("x2" in variables) and ("x3" in variables):
        pol = pol(x, y, z, 0)
    elif ("x1" in variables) and ("x2" in variables) and ("x4" in variables):
        pol = pol(x, y, 0, z)
    elif ("x1" in variables) and ("x3" in variables) and ("x4" in variables):
        pol = pol(x, 0, y, z)
    elif ("x2" in variables) and ("x3" in variables) and ("x4" in variables):
        pol = pol(0, x, y, z)
    
    xoffset = 0
    while pol(xoffset,0,0) == 0:
        xoffset += 1
    pol = pol(x + xoffset, y, z)
    while gcd(pol(0, 0, 0), X) != 1:
        X = next_prime(X, proof=False)
    while gcd(pol(0, 0, 0), Y) != 1:
        Y = next_prime(Y, proof=False)
    while gcd(pol(0, 0, 0), Z) != 1:
        Z = next_prime(Z, proof=False)

    pol = P(pol/gcd(pol.coefficients()))
    p000 = pol(0, 0, 0)

    delta = max(pol.degree(x), pol.degree(y), pol.degree(z))

    W = max(abs(i) for i in pol(x*X,y*Y,z*Z).coefficients())
    u = W + ((1-W) % abs(p000))
    N = u*(X*Y*Z)^l # modulus for polynomials

    p000inv = inverse_mod(p000,N)
    polq = P(sum((i*p000inv % N)*j for i, j in zip(pol.coefficients(), pol.monomials())))
    polynomials = []
    for i in range(delta+l+1):
        for j in range(delta+l+1):
            for k in range(delta+l+1):
                # if 0 <= i <= l and 0 <= j <= l and 0 <= k <= l:
                if 0 <= i + j + k <= l:
                    polynomials.append(polq * x^i * y^j * z^k * X^(l-i) * Y^(l-j) * Z^(l-k))
                else:
                    polynomials.append(x^i * y^j * z^k * N)

    monomials = []
    for i in polynomials:
        for j in i.monomials():
            if j not in monomials:
                monomials.append(j)
    monomials.sort()

    L = matrix(ZZ,len(monomials))
    for i in range(len(monomials)):
        for j in range(len(monomials)):
            L[i,j] = polynomials[i](X*x,Y*y,Z*z).monomial_coefficient(monomials[j])

    L = matrix(ZZ,sorted(L, reverse=True))
    L = L.LLL()

    roots = []
    P2.<q> = PolynomialRing(ZZ)

    for i in range(L.nrows()-1):
        for j in range(i+1,L.nrows()):
            pol2 = P(sum(map(mul, zip(L[i],monomials)))(x/X,y/Y,z/Z))
            pol3 = P(sum(map(mul, zip(L[j],monomials)))(x/X,y/Y,z/Z))

            r = pol.resultant(pol2, z)
            r2 = pol.resultant(pol3, z)
            r = r.resultant(r2,y)
            assert r.is_univariate()

            if r.is_constant():
                continue

            r = r(q,0,0)

            if len(r.roots()) > 0:
                for x0, _ in r.roots():
                    if x0 == 0:
                        continue
                    for y0, _ in P2(r2(x0,q,0)).roots():
                        for z0, _ in P2(pol(x0,y0,q)).roots():
                            if pol(x0-xoffset,y0,z0) == 0:
                                return [{variables[0]: x0-xoffset, variables[1]: y0, variables[2]: z0}]
