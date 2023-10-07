load('util.sage')

def coron_trivariate(pol, X, Y, Z, l=2, sol=[], debug=True):
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

    origAlg = True

    if not origAlg:
        polynomials = []
        i, j, k = 0, 0, 0
        polynomials.append(polq * x^i * y^j * z^k * X^(l-i) * Y^(l-j) * Z^(l-k))
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
                    L[i,j] = p(X*x,Y*y,Z*z).monomial_coefficient(mon)
                else:
                    L[i,j] = N
    else:
        polynomials = []
        for i in range(delta+l+1):
            for j in range(delta+l+1):
                for k in range(delta+l+1):
                    # if 0 <= i <= l and 0 <= j <= l and 0 <= k <= l:
                    if 0 <= i + j + k <= l:
                        polynomials.append(polq * x^i * y^j * z^k * X^(l-i) * Y^(l-j) * Z^(l-k))
                    # else:
                    elif i + j + k <= delta + l:
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

    if debug:
        print("is polynomial longer than guaranteed shortest: ", is_poly_longer_than_guaranteed_shortest(L, monomials, X, Y, Z, debug=debug))
    
    L = matrix(ZZ,sorted(L, reverse=True))

    if debug:
        print("")
        print("Bitlengths of matrix elements (before reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())
    
    L = L.LLL()

    if debug:
        print("")
        print("Bitlengths of matrix elements (after reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())

    roots = []
    P2.<q> = PolynomialRing(ZZ)

    pol2 = P(sum(map(mul, zip(L[0],monomials)))(x/X,y/Y,z/Z))
    pol3 = P(sum(map(mul, zip(L[1],monomials)))(x/X,y/Y,z/Z))

    if debug:
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

        c1 = check_poly_3(pol2, N, X, Y, Z)
        c2 = check_poly_3(pol3, N, X, Y, Z)
        print("Check Howgrave-Graham lemma for the first two polynomials: ", c1, c2)

    coeffs = pol.coefficients()
    coeffs2 = pol2.coefficients()

    if len(pol2.coefficients()) == 1:
        if len(pol3.coefficients()) != 1:
            coeffs2 = pol3.coefficients()
        else:
            print("TODO")
            sys.exit(0)

    # stupid comparison, but simply pol == -pol doesn't work because the constant can be moved by N
    is_same = coeffs[0] == coeffs2[0] and coeffs[1] == coeffs2[1] and coeffs[2] == coeffs2[2] and coeffs[3] == coeffs2[3] 
    is_neg = coeffs[0] == -coeffs2[0] and coeffs[1] == -coeffs2[1] and coeffs[2] == -coeffs2[2] and coeffs[3] == -coeffs2[3] 

    """
    if is_same or is_neg:
        print("??????????????? same")
        return [], "short same"
    """

    for i in range(L.nrows()-1):
        for j in range(i+1,L.nrows()):
            pol2 = P(sum(map(mul, zip(L[i],monomials)))(x/X,y/Y,z/Z))
            pol3 = P(sum(map(mul, zip(L[j],monomials)))(x/X,y/Y,z/Z))

            pol2_evals_to_0 = pol2(sol[0], sol[1], sol[2]) == 0
            pol3_evals_to_0 = pol3(sol[0], sol[1], sol[2]) == 0

            if debug:
                print("--------------")
                print(pol2)
                print("--")
                print(pol3)
                c1 = check_poly_3(pol2, N, X, Y, Z)
                c2 = check_poly_3(pol3, N, X, Y, Z)
                print("Check Howgrave-Graham lemma for the first two polynomials: ", c1, c2)
                print(pol2(sol[0], sol[1], sol[2])) # known component is removed in small_roots
                print(pol3(sol[0], sol[1], sol[2]))

            r = pol.resultant(pol2, z)
            r2 = pol.resultant(pol3, z)
            r = r.resultant(r2,y)
            assert r.is_univariate()

            if r.is_constant():
                if debug:
                    "r is constant"
                continue

            r = r(q,0,0)
            if debug:
                print(r)
                print("len roots: ", len(r.roots()))

            if len(r.roots()) > 0:
                for x0, _ in r.roots():
                    if x0 == 0:
                        continue
                    for y0, _ in P2(r2(x0,q,0)).roots():
                        for z0, _ in P2(pol(x0,y0,q)).roots():
                            if pol(x0-xoffset,y0,z0) == 0:
                                if debug:
                                    print("found: ", i, j)
                                    print("")
                                    """
                                    print(pol)
                                    print("")
                                    print(pol2)
                                    print("")
                                    print(pol3)
                                    print("")
                                    """
                                if (pol2_evals_to_0 == 0) or (pol3_evals_to_0 == False):
                                    print("Trying to see whether sometimes succeeds without eval to 0")

                                    print(pol2(sol[0], sol[1], sol[2])) # known component is removed in small_roots
                                    print(pol3(sol[0], sol[1], sol[2]))
                                    sys.exit(0)
                                return [{variables[0]: x0-xoffset, variables[1]: y0, variables[2]: z0}], ""
    return [], "dependent"
    
    
