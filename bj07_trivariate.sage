def bj07_trivariate(pol, X, Y, Z, opts={}):
    """
    Returns all small roots of pol.

    Applies the Bauer-Joux extension of Coppersmith's method to find small
    roots of a trivariate polynomial.
    Published: https://dx.doi.org/10.1007/978-3-540-72540-4_21
    Eprint: https://iacr.org/archive/eurocrypt2007/45150361/45150361.pdf

    Args:
        pol: The polynomial to find small integer roots of.
        X: Upper limit on x.
        Y: Upper limit on y.
        Z: Upper limit on z.
        opts: Various options for debugging and fine tuning the algorithm.
            - checkFirst (True): check entire basis of lattice L1' before
              moving on to Step 2 (instead of just last element)
            - checkSecond (True): check entire basis of lattice L2'
            - returnFirst (True): return right after checking L1' if any zeros
              were found (has no effect if checkFirst=False)
            - skipGB (False): instead of computing the GB for L2, append p1 to
              existing shifts and repeat step 1
            - S ("predef0"): choose the set S from predef0-2, or specify a
              number to generate all monomials up to that degree
            - solution (None): a known root
            - debug (False): print everything

    Returns:
        A list of successfully found roots [(x0, y0, z0), ...].

    Raises:
        ValueError: If pol is not bivariate
    """

    # parse tuning options
    default_opts = {"checkFirst": True,
                    "checkSecond": True,
                    "returnFirst": True,
                    "skipGB": False,
                    "S": "predef0",
                    "solution": None,
                    "debug": False,}
    options = {}
    if not opts:
        options = default_opts
    for o in default_opts:
        if o in opts:
            options[o] = opts[o]
        else:
            options[o] = default_opts[o]
    # turn into flags
    debug = options["debug"]
    check_first = options["checkFirst"]
    check_second = options["checkSecond"]
    return_first = options["returnFirst"]
    skip_GB = options["skipGB"]
    S_set = options["S"]
    sols = options["solution"]

    if pol.nvariables() != 3:
        raise ValueError("pol is not trivariate")

    P.<x,y,z> = PolynomialRing(ZZ, order="deglex")

    variables = [str(pol.variables()[0]), str(pol.variables()[1]), str(pol.variables()[2])]

    if "x1" not in variables:
      pol = pol(0, x, y, z)
    elif "x2" not in variables:
      pol = pol(x, 0, y, z)
    elif "x3" not in variables:
      pol = pol(x, y, 0, z)
    elif "x4" not in variables:
      pol = pol(x, y, z, 0)
    
    # # W is the inf norm of pol(Xx,Yy,Zz)
    # W = max(abs(i) for i in pol(x*X,y*Y,z*Z).coefficients())
    # # d_x, d_y, d_z are max degrees
    # dx, dy, dz = pol.degrees()

    roots = []

    # STEP 1: find p1 with randomization ideal S
    # monomials of pol
    m_pol = pol.monomials()

    # the set S -- TODO: fuzz
    S = []
    if S_set == "predef0":
        S = [P(1), x*y, x*z, y*z]
    elif S_set == "predef1":
        S = [P(1), x, y, z]
    elif S_set == "predef2":
        S = [P(1), x, y, z, x*y, x*z, y*z]
    elif type(S_set) == Integer:
        l = S_set
        for i in range(0, l+1):
            for j in range(0, l+1):
                for k in range(0, l+1):
                    if i + j + k <= l:
                        S.append(x^i * y^j * z^k)
    else:
        raise ValueError("S must be either predef[0-2] or an integer")
    if debug:
        print("S: ", S)

    # constuct M from above s.t. (S, M) is admissible for pol
    M = [m1*m2 for m1 in m_pol for m2 in S]

    # # bounds s_x, s_y, s_z
    # sx = 0
    # sy = 0
    # sz = 0
    # for mon in M:
    #     if mon not in S:
    #         ds = mon.degrees()
    #         sx += ds[0]
    #         sy += ds[1]
    #         sz += ds[2]

    # add monomials of pol into M, in case 1 not in S
    for mon in m_pol:
        M.append(mon)
    M.append(P(1))
    # remove duplicates and sort
    M = list(set(M))
    M.sort()
    # define length parameters
    m, s = len(M), len(S)

    # # find the const c
    # c = ceil((m-s)^2 / (s * (dx^2 + dy^2 + dz^2)))
    # if debug:
    #     print("W = ", W)
    #     print("c = ", c)
    #     print("X, Y, Z = ", X, Y, Z)
    #     print("sx, sy, sz = ", sx, sy, sz)
    #     print("dx, dy, dz = ", dx, dy, dz)
    #     print("m, s = ", m, s)
    #     print("LHS = ", round(X^sx * Y^sy * Z^sy))
    #     print("RHS = ", round(W^s * 2^(-(6+c)*s*(dx^2 + dy^2 + dz^2))))
    # # this isn't satisfied
    # # assert X^sx * Y^sy * Z^sy < W^s * 2^(-(6+c)*s*(dx^2 + dy^2 + dz^2))

    # construct polynomials
    polys = [pol * p for p in S]

    # construct lattice L1
    # I1 has inverse bounds on the diagonal
    I1 = matrix(QQ, m)
    for i in range(m):
        ds = M[i].degrees()
        I1[i,i] = X^(-ds[0]) * Y^(-ds[1]) * Z^(-ds[2])
    # P1 has S*pol as columns
    # WARNING, has to be ZZ, not QQ, otherwise Sage has a very liberal
    # interpretation of what a unimodular transformation means
    P1 = matrix(ZZ, m, s)
    for j in range(s):
        for i in range(m):
            P1[i, j] = polys[j].monomial_coefficient(M[i])

    # find sublattice L1' of L1 with all zeroes in last s columns
    P1p, U1 = P1.hermite_form(transformation=True)
    # assert P1p == U1 * P1
    I1p = U1 * I1
    # assert I1p != I1
    if debug:
        set_verbose(2)
        print("P1 matrix elements:")
        print(P1)
        print("P1p matrix elements:")
        print(P1p)
        print("U1 matrix elements:")
        print(U1)
    if debug:
        print("I1p matrix elements (before reduction):")
        print(I1p)
    # extract bottom m-s elements
    L1pr = I1p[s:]

    # L1' ---LLL--> B ---GS--> B*
    L1pr = L1pr.LLL()
    if debug:
        print("I1p matrix elements (after reduction):")
        print(L1pr)
    GS1, _ = L1pr.gram_schmidt()
    if debug:
        print("I1p matrix elements (after GS):")
        print(GS1)
    # check solutions -- DELETE
    if debug and sols is not None:
        print("=== BEGIN STEP 1 polynomials on solution")
        print(pol(*sols))
        PP.<xx, yy, zz> = PolynomialRing(QQ)
        MM = []
        for mon in M:
            ds = mon.degrees()
            MM.append(xx^ds[0] * yy^ds[1] * zz^ds[2])
        x0, y0, z0 = sols
        s = []
        for mon in MM:
            ds = mon.degrees()
            s.append((x0 / X)^ds[0] * (y0 / Y)^ds[1] * (z0 / Z)^ds[2])
        s0 = vector(QQ, s)
        s0norm2 = s0.norm().numerical_approx()
        s0norm1 = s0.norm(1).numerical_approx()
        s0normInf = s0.norm(Infinity).numerical_approx()
        print("NORM2 s0 =", s0norm2)
        print("NORM1 s0 =", s0norm1)
        print("NORMInf s0 =", s0normInf)
        for b in GS1:
            bnorm2 = b.norm().numerical_approx()
            bnorm1 = b.norm(1).numerical_approx()
            bnormInf = b.norm(Infinity).numerical_approx()
            print("NORM2 bi =", bnorm2, bnorm2 >= s0norm2)
            print("NORM1 bi =", bnorm1, bnorm1 >= s0norm1)
            print("NORMInf bi =", bnormInf, bnormInf >= s0normInf)
            p = PP(sum(map(mul, zip(b,MM))))
            print(p(sols[0]/X, sols[1]/Y, sols[2]/Z), p)
            print(p(xx/X, yy/Y, zz/Z)(sols[0], sols[1], sols[2]))
            p = p(xx/X, yy/Y, zz/Z)
            p_lcm = lcm(list(map(lambda x: x.denominator(), p.coefficients())))
            p = p_lcm * p
            p = P(p)
            print(p(*sols), p)
        print("=== END STEP 1 polynomials on solution")

    # extract b_r* as an integer polynomial (see above warning on
    # unimodularity), i.e. we scale it by LCM of coefficient denominators,
    # which preserves roots
    if debug:
        print("br* = ", GS1[-1])
    coefs1 = []
    for i in range(m):
        coefs1.append(GS1[-1][i]/M[i](X,Y,Z))
    lcm1 = lcm([coef.denominator() for coef in coefs1])
    coefs1 = [lcm1 * coef for coef in coefs1]
    p1 = P(sum(map(mul, zip(coefs1,M))))
    if debug and sols is not None:
        print("p1: ", p1(*sols), p1)
    if debug:
        print("STEP 1 -- p1: ", p1)

    # IDEA: Try b_{r-1}* as well, maybe we got lucky and produced 2+
    # polynomials in the first step if the solution is extra short. In fact why
    # not just for loop over them in reverse order.
    # try:
    if check_first:
        for j in range(GS1.nrows()-2, -1, -1):
            if debug:
                # print("trying b[r-1]* = ", GS1[-2])
                print("trying b[{}]* = ".format(j), GS1[j])
            coefs1p = []
            for i in range(m):
                # coefs1p.append(GS1[-2][i]/M[i](X,Y,Z))
                coefs1p.append(GS1[j][i]/M[i](X,Y,Z))
            lcm1p = lcm([coef.denominator() for coef in coefs1p])
            coefs1p = [lcm1p * coef for coef in coefs1p]
            p1p = P(sum(map(mul, zip(coefs1p,M))))
            if debug and sols is not None:
                print("p1p: ", p1p(*sols), p1p)
            r1 = pol.resultant(p1, z)
            r2 = pol.resultant(p1p, z)
            r3 = r1.resultant(r2, y)
            if debug:
                print("STEP 1.5 -- r1: ", r1)
                print("STEP 1.5 -- r2: ", r2)
                print("STEP 1.5 -- r3: ", r3)
            # extract univariate polynomial zeros
            if r3.is_constant():
                if debug:
                    print("STEP 1.5 -- r3 is constant!")
                continue
            P2.<q> = PolynomialRing(ZZ)
            r3 = r3(q, 0, 0)
            if len(r3.roots()) == 0:
                if debug:
                    print("STEP 1.5 -- r3 has no roots!")
                continue
            for x0, _ in r3.roots():
                if debug:
                    print("Potential x0: ", x0)
                r3y1 = P2(r1(x0, q, 0))
                if r3y1.is_constant():
                    continue
                for y0, _ in r3y1.roots():
                    if debug:
                        print("Potential y0: ", y0)
                    r3z = P2(pol(x0, y0, q))
                    if r3z.is_constant():
                        continue
                    for z0, _ in r3z.roots():
                        if debug:
                            print("Potential z0: ", z0)
                            print("Found root: ({}, {}, {})\n".format(x0, y0, z0))
                        if {variables[0]: x0, variables[1]: y0, variables[2]: z0} not in roots:
                            roots.append({variables[0]: x0, variables[1]: y0, variables[2]: z0})
                r3y2 = P2(r2(x0, q, 0))
                if r3y2.is_constant():
                    continue
                for y0, _ in r3y2.roots():
                    if debug:
                        print("Potential y0: ", y0)
                    r3z = P2(pol(x0, y0, q))
                    if r3z.is_constant():
                        continue
                    for z0, _ in r3z.roots():
                        if debug:
                            print("Potential z0: ", z0)
                            print("Found root: ({}, {}, {})\n".format(x0, y0, z0))
                        if {variables[0]: x0, variables[1]: y0, variables[2]: z0} not in roots:
                            roots.append({variables[0]: x0, variables[1]: y0, variables[2]: z0})
                        # if we have a root just return?
                        if return_first and len(roots) > 0:
                            return roots
        if debug and len(roots) == 0:
            print("First step did not yield extra polynomials.")
        # if we have a root just return?
        if return_first and len(roots) > 0:
            return roots

    # STEP 2: find p2 with truncated GB of I = (pol, p1)
    # auxiliary function to check p defined over M
    def pol_defined_over(poly, mons):
        for m in poly.monomials():
            if m not in mons:
                return False
        return True

    I = ideal(pol, p1)

    F = []
    if skip_GB:
        # use heuristic F = [S*pol, p1]
        F = [p for p in polys]
        F.append(p1)
    else:
        # compute truncated GB
        G = list(I.groebner_basis())
        if debug:
            print("STEP 2 -- GB: ", len(G), G)
        trG = [g for g in G if pol_defined_over(g, M)]
        if debug:
            print("trGB: ", len(trG), trG)
        assert len(trG) > 0
        # all non-const monomials up to n_bound -- potentially unnecessary
        n_bound = M[-1].degree()
        M_n = []
        for i in range(0, n_bound):
            for j in range(0, n_bound):
                for k in range(0, n_bound):
                    if i + j + k <= n_bound:
                        M_n.append(x^i * y^j * z^k)
        # add multiples of GB elements up to degree bound -- potentially unneccessary
        mtrG = [g * mn for g in trG for mn in M if pol_defined_over(g*mn, M)]
        if debug:
            print("trGB after mult: ", len(mtrG), mtrG)
        F = mtrG
        if len(mtrG) >= m:
            # use heuristic F = [S*pol, p1] anyway
            F = [p for p in polys]
            F.append(p1)
    assert len(mtrG) != 0
    # assert len(mtrG) < m
    
    # define length parameters
    # m = len(M)
    t = len(mtrG)

    # new lattice L2
    # I2 is the same as I1, inverse bounds on diagonal
    I2 = matrix(QQ, I1)
    # P2 has multiplied truncated GB as columns (again must be ZZ, see P1
    # warning)
    P2 = matrix(ZZ, m, t)
    for j in range(t):
        for i in range(m):
            P2[i, j] = mtrG[j].monomial_coefficient(M[i])
    P2p, U2 = P2.hermite_form(transformation=True)
    assert P2p == U2 * P2
    I2p = U2*I2
    assert I2p != I2
    if debug:
        set_verbose(2)
        print("P2 matrix elements:")
        print(P2)
        print("P2p matrix elements:")
        print(P2p)
        print("U2 matrix elements:")
        print(U2)
    if debug:
        print("I2p matrix elements (before reduction):")
        print(I2p)
    # extract bottom m-t elements
    L2pr = I2p[t:]

    # L2' ---LLL--> B ---GS--> B*
    L2pr = L2pr.LLL()
    if debug:
        print("I2p matrix elements (after reduction):")
        print(L2pr)
    GS2, _ = L2pr.gram_schmidt()
    if debug:
        print("I2p matrix elements (after GS):")
        print(GS2)
    # check solutions -- DELETE
    if debug and sols is not None:
        print("=== BEGIN STEP 2 polynomials on solution")
        print(pol(*sols), p1(*sols))
        PP.<xx, yy, zz> = PolynomialRing(QQ)
        MM = []
        for mon in M:
            ds = mon.degrees()
            MM.append(xx^ds[0] * yy^ds[1] * zz^ds[2])
        for b in GS2:
            p = PP(sum(map(mul, zip(b,MM))))
            print(p(sols[0]/X, sols[1]/Y, sols[2]/Z), p)
            print(p(xx/X, yy/Y, zz/Z)(sols[0], sols[1], sols[2]))
            p = p(xx/X, yy/Y, zz/Z)
            p_lcm = lcm(list(map(lambda x: x.denominator(), p.coefficients())))
            p = p_lcm * p
            p = P(p)
            print(p(*sols), p)
        print("=== END STEP 2 polynomials on solution")

    # extract c_r* as an integer polynomial (see above warning on
    # unimodularity), i.e. we scale it by LCM of coefficient denominators,
    # which preserves roots
    lbound2 = GS2.nrows() - 2
    if check_second:
        lbound2 = -1
    for j in range(GS2.nrows()-1, lbound2, -1):
        if debug:
            print("c[{}]* = ".format(j), GS2[j])
        coefs2 = []
        for i in range(m):
            coefs2.append(GS2[j][i]/M[i](X,Y,Z))
        lcm2 = lcm([coef.denominator() for coef in coefs2])
        coefs2 = [lcm2 * coef for coef in coefs2]
        p2 = P(sum(map(mul, zip(coefs2,M))))
        if debug and sols is not None:
            print("p2: ", p2(*sols), p2)
        if debug:
            print("STEP 2 -- p2: ", p2)
        # extract polynomials
        P2.<q> = PolynomialRing(ZZ)
        if debug:
            print("STEP 2 -- p2: ", p2)
            print("STEP 2 -- p2 in I?: ", p2 in I)

        # STEP 3 -- find solutions
        # Macaulay resultant seems to only be supported for homogeneous
        # polynomials, which would mean introducing another variable...
        # compute resultant tree instead
        r1 = pol.resultant(p1, z)
        r2 = pol.resultant(p2, z)
        r3 = r1.resultant(r2, y)
        if debug:
            print("STEP 3 -- r1: ", r1)
            print("STEP 3 -- r2: ", r2)
            print("STEP 3 -- r3: ", r3)
        # extract univariate polynomial zeros
        if r3.is_constant():
            if debug:
                print("STEP3 -- r3 is constant!")
            continue
        r3 = r3(q, 0, 0)
        if len(r3.roots()) == 0:
            if debug:
                print("STEP 3 -- r3 has no roots!")
            continue
        for x0, _ in r3.roots():
            if debug:
                print("Potential x0: ", x0)
            r3y1 = P2(r1(x0, q, 0))
            if r3y1.is_constant():
                continue
            for y0, _ in r3y1.roots():
                if debug:
                    print("Potential y0: ", y0)
                r3z = P2(pol(x0, y0, q))
                if r3z.is_constant():
                    continue
                for z0, _ in r3z.roots():
                    if debug:
                        print("Potential z0: ", z0)
                        print("Found root: ({}, {}, {})\n".format(x0, y0, z0))
                    if {variables[0]: x0, variables[1]: y0, variables[2]: z0} not in roots:
                        roots.append({variables[0]: x0, variables[1]: y0, variables[2]: z0})
            r3y2 = P2(r2(x0, q, 0))
            if r3y2.is_constant():
                continue
            for y0, _ in r3y2.roots():
                if debug:
                    print("Potential y0: ", y0)
                r3z = P2(pol(x0, y0, q))
                if r3z.is_constant():
                    continue
                for z0, _ in r3z.roots():
                    if debug:
                        print("Potential z0: ", z0)
                        print("Found root: ({}, {}, {})\n".format(x0, y0, z0))
                    if {variables[0]: x0, variables[1]: y0, variables[2]: z0} not in roots:
                        roots.append({variables[0]: x0, variables[1]: y0, variables[2]: z0})
    return roots

# if __name__ == "__main__":
#     P.<x1, x2, x3, x4> = PolynomialRing(ZZ)
#
#     p = 11956566944641502957704189594909498993478297403838643406058180334130656750161
#     Q.<i,j,k> = QuaternionAlgebra(QQ, -1, -p)
#     I = Q.ideal([2, 1 - j, i + j, 3/2 + 1/2*i + 1/2*j + 1/2*k])
#
#     iso_degree = 2**160
#     quadratic_form = 120864756928231872170197654963601902031*x1^2 + 40379784023742363682691476922140013988*x1*x2 - 80210377761494289609629402238643748110*x1*x3 - 49019009704082514622267696968581725297*x1*x4 + 194680861828862664977080049082731556942*x2^2 + 156088821196201368396536230067211356461*x2*x3 + 147083019515271147858011236632623029956*x2*x4 + 410215661800567766037607854071586163516*x3^2 + 205107830900283883018803927035793081758*x3*x4 + 615323492700851649056411781107379245274*x4^2
#     sol = [99533, -30399, 42217, -4004]
#     known_variables = [2]
#     sol_x = [99533, -30399, -4004]
#
#     b = 5000000000000
#     b = 1661926 # max for the above example and S = [1, x, y, z, x*y, x*z, y*z]
#     b = 1693442 # max for the above example and S = [1, x, y, z]
#     b = 1258884 # max for the above example and S = [1, x*y, y*z, x*z]
#     b = 14956000 # max for S deg up to 2
#     b = 10000000
#
#     f = quadratic_form - iso_degree 
#     for v in known_variables:
#       if v == 0:
#         f = f(x1 = sol[0])
#       elif v == 1:
#         f = f(x2 = sol[1])
#       elif v == 2:
#         f = f(x3 = sol[2])
#       elif v == 3:
#         f = f(x4 = sol[3])
#
#     roots = bj07_trivariate(f, b, b, b, opts={"debug":True, "S": 2, "solution": sol_x})
#     print(roots)
