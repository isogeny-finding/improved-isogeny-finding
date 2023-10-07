"""
p is a polynomial in
P.<x1, x2, x3, x4> = PolynomialRing(ZZ)
p = a * x^2 + b * x * y + c * y^2 + d * x + e * y + f = 0

Uses method from:
A New Look at an Old Equation
by R.E. Sawilla, A.K. Silvester, and H.C. Williams
"""
def transform_to_cornacchia_and_solve(p):
  var1 = p.variables()[0]
  var2 = p.variables()[1]

  a = p.monomial_coefficient(var1^2)
  b = p.monomial_coefficient(var1 * var2)
  c = p.monomial_coefficient(var2^2)
  d = p.monomial_coefficient(var1)
  e = p.monomial_coefficient(var2)
  f = p(0, 0, 0, 0)

  D = b^2 - 4 * a * c
  E = b * d - 2 * a * e
  F = d^2 - 4 * a * f
  N = E^2 - D * F

  # Is D always < 0?
  assert D < 0

  # D * Y^2 = (D * y + E)^2 + D * F - E^2
  # X = D * y + E
  # Y = 2 * a * x + b * y + d
  Y = 2 * a * var1 + b * var2 + d
 
  print("D: ", D)
  print("E: ", E)
  print("N: ", N)

  Q = BinaryQF([1, 0, -D])
  sol = Q.solve_integer(N)
  print(sol)

  y = (sol[0] - E) / D
  x = (sol[1] - b * y - d) / (2 * a)
  print(x, y)
  check = p.subs({var1: x, var2: y}) # p(var1 = x) doesn't work
  assert check == 0
  return x, y

"""
P.<x1, x2, x3, x4> = PolynomialRing(ZZ)
p = 4*x1^2 + 4*x2^2 + 6*x2
transform_to_cornacchia_and_solve(p)
"""