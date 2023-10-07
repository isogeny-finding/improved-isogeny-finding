def norm2_square(p):
    if type(p) == Integer:
        return p*p
    n = 0
    for c in p.coefficients():
        n += c*c
    return n

def is_poly_longer_than_guaranteed_shortest(L, monomials, X, Y, Z=None, T=None, debug=True):
    """
    L is the matrix where rows represent the vectors.
    The LLL ensures that the shortest vector it finds is smaller than `bound = 2**((nn-1)/4) * dt**(1/nn)`.
    We check whether the first row of L is bigger than `bound`.
    If it is, the LLL will return find a smaller vector.
    """

    if Z == None and T == None:
        P.<x,y> = PolynomialRing(ZZ)
        pp = P(sum(map(mul, zip(L[0],monomials)))(x,y))
    elif T == None:
        P.<x,y,z> = PolynomialRing(ZZ)
        pp = P(sum(map(mul, zip(L[0],monomials)))(x,y,z))
    else:
        P.<x,y,z,t> = PolynomialRing(ZZ)
        pp = P(sum(map(mul, zip(L[0],monomials)))(x,y,z,t))
    dt = L.determinant().abs()
    nn = L.dimensions()[0]
    # ∥b1∥ ≤ 2^((ω−1)/4) det(L)^(1/ω)
    n2 = norm2_square(pp)

    # bound^(4 * nn) = 2**((nn-1) * nn) * dt**4
    bound_exp = 2**((nn-1) * nn) * dt**4
    n2_exp = n2**(2 * nn)

    print("")
    print(n2_exp.bit_length())
    print(bound_exp.bit_length())

    # we want n2_exp > bound_exp, because then we have a guarantee that LLL will shorten it 
    return n2_exp > bound_exp

def check_poly_2(p, N, X, Y):
    # omega * norm(p(xX, yY))**2 < n**2 

    P.<x,y,z,t> = PolynomialRing(ZZ)
    omega = len(p.coefficients())
    norm_square = norm2_square(p(X*x, Y*y))

    return omega * norm_square < N*N

def check_independence_with_poly_2(p, W, k, X, Y):
    P.<x,y,z,t> = PolynomialRing(ZZ)
    norm_square = norm2_square(p(X*x, Y*y))
    # ∥h(xX,yY)∥<2^(−ω) ·(XY)^k ·W
    # 2^-(d+1)^2
    # rhs = 2^(-9) * (X * Y)**k * W
    d = 2
    rhs = 2^(-(k + d + 1)*(k + d + 2) / 2) * (X * Y)**k * W
    # rhs = 2^(-(k + d + 1)**2) * (X * Y)**k * W

    return norm_square < rhs**2

def check_independence_with_det_2(L, W, k, X, Y):
    P.<x,y,z,t> = PolynomialRing(ZZ)
    # 2^((ω−1)/4) det(L)^(1/ω) <2^(−ω) ·(XY)^k ·W
    # rhs = 2^(-9) * (X * Y)**k * W
    d = 2
    omega = (k + d + 1)*(k + d + 2) / 2
    # omega = (k + d + 1)*(k + d + 1)

    dt = L.determinant().abs()
    nn = L.dimensions()[0]
    lhs = 2**((nn-1) * nn) * dt**4

    assert omega == nn

    rhs1 = 2**(-omega) * (X * Y)**k * W
    rhs = rhs1**(4*nn)

    return lhs < rhs

def check_poly_3(p, N, X, Y, Z):
    # omega * norm(p(xX, yY, zZ))**2 < n**2 

    P.<x,y,z,t> = PolynomialRing(ZZ)
    omega = len(p.coefficients())
    norm_square = norm2_square(p(X*x, Y*y, Z*z))

    return omega * norm_square < N*N

def check_poly_4(p, N, X, Y, Z, T):
    # omega * norm(p(xX, yY, zZ, tT))**2 < n**2 

    P.<x,y,z,t> = PolynomialRing(ZZ)
    omega = len(p.coefficients())
    norm_square = norm2_square(p(X*x, Y*y, Z*z, T*t))

    return omega * norm_square < N*N
