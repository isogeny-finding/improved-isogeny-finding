def coron(pol, X, Y, k=2, sol=[], debug=True):
    """
    Returns the first root of pol that is found.

    Applies the Coron's generalization of Coppersmith's algorithm for finding small
    integer roots of bivariate polynomials modulo an integer.

    Coron, J.S.: Finding small roots of bivariate integer polynomial equations:
    A direct approach. In: CRYPTO. vol. 4622, pp. 379–394. Springer (2007)

    The implementation slightly adapted from https://github.com/ubuntor/coppersmith-algorithm.

    Args:
        pol: The polynomial to find small integer roots of.
        X: Upper limit on x.
        Y: Upper limit on y.
        k: Determines size of lattice.

    Returns:
        The first root that is found.

    Raises:
        ValueError: If pol is not bivariate
    """

    if pol.nvariables() != 2:
        raise ValueError("pol is not bivariate")

    P.<x,y> = PolynomialRing(ZZ)

    variables = [str(pol.variables()[0]), str(pol.variables()[1])]
    if ("x1" in variables) and ("x2" in variables):
        pol = pol(x, y, 0, 0)
    elif ("x1" in variables) and ("x3" in variables):
        pol = pol(x, 0, y, 0)
    elif ("x1" in variables) and ("x4" in variables):
        pol = pol(x, 0, 0, y)
    elif ("x2" in variables) and ("x3" in variables):
        pol = pol(0, x, y, 0)
    elif ("x2" in variables) and ("x4" in variables):
        pol = pol(0, x, 0, y)
    elif ("x3" in variables) and ("x4" in variables):
        pol = pol(0, 0, x, y)

    xoffset = 0
    while pol(xoffset, 0) == 0:
        xoffset += 1

    pol = pol(x + xoffset, y)
    while gcd(pol(0, 0), X) != 1:
        X = next_prime(X, proof=False)
    while gcd(pol(0, 0), Y) != 1:
        Y = next_prime(Y, proof=False)

    pol = P(pol/gcd(pol.coefficients()))
    p00 = pol(0, 0)
    delta = max(pol.degree(x), pol.degree(y))
    W = max(abs(i) for i in pol(x*X, y*Y).coefficients())

    u = W + ((1-W) % abs(p00))
    N = u*(X*Y)^k

    p00inv = inverse_mod(p00,N)
    polq = P(sum((i*p00inv % N)*j for i, j in zip(pol.coefficients(), pol.monomials())))
    polynomials = []
    for i in range(delta+k+1):
        for j in range(delta+k+1):
            # if 0 <= i <= k and 0 <= j <= k:
            if 0 <= i + j <= k:
                polynomials.append(polq * x^i * y^j * X^(k-i) * Y^(k-j))
            # else:
            elif i + j <= delta + k:
                polynomials.append(x^i * y^j * N)

    monomials = []
    for i in polynomials:
        for j in i.monomials():
            if j not in monomials:
                monomials.append(j)
    monomials.sort()

    L = matrix(ZZ,len(monomials))
    for i in range(len(monomials)):
        for j in range(len(monomials)):
            L[i, j] = polynomials[i](X*x, Y*y).monomial_coefficient(monomials[j])

    if debug:
        print("is polynomial longer than guaranteed shortest: ", is_poly_longer_than_guaranteed_shortest(L, monomials, X, Y, debug=debug))

    L = matrix(ZZ,sorted(L, reverse=True)) 

    if debug:
        print("")
        print("Bitlengths of matrix elements (before reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())

        """
        pp = P(sum(map(mul, zip(L[0],monomials)))(x,y))
        orig_norm = norm2_square(pp)
        print("")
        print(pp)
        print("")
        print("orig norm: ")
        print(orig_norm)
        print("")
        """

    L = L.LLL(delta=0.75)
    roots = []

    pol2 = P(sum(map(mul, zip(L[0],monomials)))(x/X,y/Y))
    if debug:
        print("Bitlengths of matrix elements (after reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())

        """
        pp = P(sum(map(mul, zip(L[0],monomials)))(x,y))
        shortened_norm = norm2_square(pp)
        print("")
        print(pp)
        print("")
        print("shortened norm: ")
        print(shortened_norm)
        print("")
        """

        """
        print("")
        print("original polynomial")
        print(pol)
        print("")
        print("first found polynomial")
        print(pol2)
        print("")

        c1 = check_poly_2(pol2, N, X, Y)
        print("Check Howgrave-Graham lemma for the first polynomial: ", c1)
        """

    coeffs = pol.coefficients()
    coeffs2 = pol2.coefficients()
    # stupid comparison, but simply pol == -pol doesn't work because the constant can be moved by N
    is_same = coeffs[0] == coeffs2[0] and coeffs[1] == coeffs2[1] and coeffs[2] == coeffs2[2] and coeffs[3] == coeffs2[3] 
    is_neg = coeffs[0] == -coeffs2[0] and coeffs[1] == -coeffs2[1] and coeffs[2] == -coeffs2[2] and coeffs[3] == -coeffs2[3] 

    """
    if is_same or is_neg:
        if debug:
            print("same")
        return [], "short same"
    """

    if debug:
        print("original polynomial")
        print(polq)
        print("")
        print(polq(x*X, y*Y))
        print("")
        c = X * Y * polq(x*X, y*Y)
        print(c)
        print("")
        l = map(lambda x: x.nbits(), c.coefficients())
        print(list(l)[::-1])
        print("")

    for i in range(L.nrows()):
        pol2 = P(sum(map(mul, zip(L[i], monomials)))(x/X, y/Y))

        if len(pol2.coefficients()) == 1:
            continue

        r = pol.resultant(pol2, y)

        if debug:
            print("pol2 (", i, ")")
            print(pol2)
            print(sol)
            print("")
            
            print(pol2 / (X*Y))
            print("")

            pol2o = P(sum(map(mul, zip(L[i], monomials)))(x, y))
            l = map(lambda x: x.nbits(), pol2o.coefficients())
            print(list(l)[::-1])
            print("")

            eval = pol2(sol[0], sol[1])
            print("eval: ", eval) # known component is removed in small_roots

        if r.is_constant():
            if debug:
                print("r is constant: ", r)
            continue

        for x0, _ in r.univariate_polynomial().roots():
            if debug:
                print("111111111111")
            if x0-xoffset in [i[0] for i in roots]:
                continue
            if debug:
                print("2222222222")
            for y0, _ in pol(x0, y).univariate_polynomial().roots():
                if debug:
                    print("3333333333333")
                if (x0-xoffset, y0) not in roots and pol(x0, y0) == 0:
                    if debug:
                        print("found: ", i)
                        print("")
                        print(pol)
                        print("")
                        print(pol2)
                        print("")

                        c1 = check_independence_with_poly_2(pol2, W, k, X, Y)
                        c2 = check_independence_with_det_2(L, W, k, X, Y)
                        print("Independence 1: ", c1)
                        print("Independence 2: ", c2)

                    roots.append({variables[0]: x0-xoffset, variables[1]: y0})
                    return roots, ""

    return [], "dependent"
