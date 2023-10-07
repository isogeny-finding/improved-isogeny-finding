def coron_quadvariate(pol, X, Y, Z, T, l=2, debug=True):
    P.<x,y,z,t> = PolynomialRing(ZZ)

    variables = [str(pol.variables()[0]), str(pol.variables()[1]), str(pol.variables()[2]), str(pol.variables()[3])]

    xoffset = 0
    while pol(xoffset, 0, 0, 0) == 0:
        xoffset += 1
    pol = pol(x+xoffset, y, z, t)
    while gcd(pol(0, 0, 0, 0), X) != 1:
        X = next_prime(X, proof=False)
    while gcd(pol(0, 0, 0, 0), Y) != 1:
        Y = next_prime(Y, proof=False)
    while gcd(pol(0, 0, 0, 0), Z) != 1:
        Z = next_prime(Z, proof=False)
    while gcd(pol(0, 0, 0, 0), T) != 1:
        T = next_prime(T, proof=False)

    pol = P(pol/gcd(pol.coefficients()))
    p0000 = pol(0,0,0,0)

    delta = max(pol.degree(x),pol.degree(y),pol.degree(z),pol.degree(t))
    W = max(abs(i) for i in pol(x*X,y*Y,z*Z,t*T).coefficients())
    u = W + ((1-W) % abs(p0000))

    if gcd(p0000, u) != 1:
        print("gcd not right")
        return []

    N = u*(X*Y*Z*T)^l
    p0000inv = inverse_mod(p0000,N)
    polq = P(sum((i*p0000inv % N)*j for i,j in zip(pol.coefficients(), pol.monomials())))

    polynomials = []
    i, j, k, u = 0, 0, 0, 0
    polynomials.append(polq * x^i * y^j * z^k * t^u * X^(l-i) * Y^(l-j) * Z^(l-k) * T^(l-u))
    # More shifted polynomials can be appended.

    mons = []
    for po in polynomials:
        mons = mons + po.monomials()

    first_mons = []
    for po in polynomials:
        ms = po.monomials()
        ms.sort()
        first_mons.append(ms[0])

    polys = []
    for i in polynomials:
        for j in i.monomials():
            if (j in first_mons) or j * N in polys:
                continue
            polys.append(j * N)

    polynomials = polynomials + polys
        
    monomials = []
    for i in polynomials:
        if type(i) == Integer:
            monomials.append(i)
            continue
        for j in i.monomials():
            if j not in monomials:
                monomials.append(j)

    monomials.sort()

    L = matrix(ZZ,len(monomials))
    for i in range(len(monomials)):
        for j in range(len(monomials)):
            mon = monomials[j]
            if type(mon) != Integer:
                mon = mon.change_ring(ZZ)
            p = polynomials[i]
            if type(p) != Integer:
                p = p.change_ring(ZZ)
            if type(mon) != Integer:
                L[i,j] = p(X*x,Y*y,Z*z,T*t).monomial_coefficient(mon)
            else:
                L[i,j] = N

    if debug:
        print("is polynomial longer than guaranteed shortest: ", is_poly_longer_than_guaranteed_shortest(L, monomials, X, Y, Z, T, debug=debug))

    L = matrix(ZZ, sorted(L, reverse=True))

    if debug:
        print("")
        print("Bitlengths of matrix elements (before reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())

        pp = P(sum(map(mul, zip(L[0],monomials)))(x,y,z,t))
        orig_norm = norm2_square(pp)
        print("")
        print(pp)
        print("")
        print("orig norm: ")
        print(orig_norm)
        print("")

    L = L.LLL()

    roots = []
    P2.<q> = PolynomialRing(ZZ)

    pol2 = P(sum(map(mul, zip(L[0],monomials)))(x/X,y/Y,z/Z,t/T))
    pol3 = P(sum(map(mul, zip(L[1],monomials)))(x/X,y/Y,z/Z,t/T))
    pol4 = P(sum(map(mul, zip(L[2],monomials)))(x/X,y/Y,z/Z,t/T))

    if debug:
        print("Bitlengths of matrix elements (after reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())

        pp = P(sum(map(mul, zip(L[0],monomials)))(x,y,z,t))
        shortened_norm = norm2_square(pp)
        print("")
        print(pp)
        print("")
        print("shortened norm: ")
        print(shortened_norm)
        print("")
    
        print("original polynomial:")
        print(pol)
        print("")
        print("first found polynomial:")
        print(pol2)
        print("")
        print("second found polynomial:")
        print(pol3)
        print("")
        print("third found polynomial:")
        print(pol4)
        print("")

        c1 = check_poly_4(pol2, N, X, Y, Z, T)
        c2 = check_poly_4(pol3, N, X, Y, Z, T)
        c3 = check_poly_4(pol4, N, X, Y, Z, T)
        print("Check Howgrave-Graham lemma for the first two polynomials: ", c1, c2, c3)

    coeffs = pol.coefficients()
    coeffs2 = pol2.coefficients()
    # stupid comparison, but simply pol == -pol doesn't work because the constant can be moved by N
    is_same = coeffs[0] == coeffs2[0] and coeffs[1] == coeffs2[1] and coeffs[2] == coeffs2[2] and coeffs[3] == coeffs2[3] 
    is_neg = coeffs[0] == -coeffs2[0] and coeffs[1] == -coeffs2[1] and coeffs[2] == -coeffs2[2] and coeffs[3] == -coeffs2[3] 

    """
    if is_same or is_neg:
        print("??????????????? same")
        return [], "short same"
    """

    for i in range(L.nrows()-1):
        pol2 = P(sum(map(mul, zip(L[i],monomials)))(x/X, y/Y, z/Z, t/T))

        if (pol == pol2) or len(pol2.monomials()) == 1:
            continue
        if pol2.nvariables() == 0:
            continue

        for j in range(i+1,L.nrows()):
            pol3 = P(sum(map(mul, zip(L[j], monomials)))(x/X, y/Y, z/Z, t/T))

            if (pol3.nvariables() == 0) or len(pol3.monomials()) == 1:
                continue

            for u in range(j+1,L.nrows()):
                pol4 = P(sum(map(mul, zip(L[u],monomials)))(x/X,y/Y,z/Z,t/T))

                if (pol4.nvariables() == 0) or len(pol4.monomials()) == 1:
                    continue
                 
                r = pol.resultant(pol2, t)
                r2 = pol.resultant(pol3, t)
                r3 = pol.resultant(pol4, t)

                rr = r.resultant(r2, z)
                rr1 = r.resultant(r3, z)
                rrr = rr.resultant(rr1, y)

                r = rrr(q,0,0,0)

                if rrr.nvariables() == 1:
                    if len(r.roots()) > 0:
                        for x0, _ in r.roots():
                            if x0 == 0:
                                continue
                            for y0, _ in P2(rr1(x0,q,0,0)).roots():
                                for z0, _ in P2(r3(x0,y0,q,0)).roots():
                                    for t0, _ in P2(pol(x0,y0,z0,q)).roots():
                                        if pol(x0-xoffset,y0,z0,t0) == 0:
                                            if debug:
                                                print("")
                                                print("found: ", i, j, k)
                                                print("")
                                                print(pol)
                                                print("")
                                                print(pol2)
                                                print("")
                                                print(pol3)
                                                print("")
                                                print(pol4)
                                                print("")
                                            return [{variables[0]: x0-xoffset, variables[1]: y0, variables[2]: z0, variables[3]: t0}], ""
                    else:
                        continue
    return [], "dependent"