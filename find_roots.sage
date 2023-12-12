import itertools
load('coron_bivariate_ana.sage')
load('coron_trivariate_ana.sage')
load('coron_quadvariate.sage')
# load('boneh_durfee_trivariate.sage')
load('bj07_trivariate.sage')

def small_roots(f, constant_polys, bounds, p, sol, m=1, d=None, method="coron", debug=True):
  if not d:
    d = f.degree()
  N = p
  assert f(x1 = sol[0], x2 = sol[1], x3 = sol[2], x4 = sol[3]) == 0
 
  if f.nvariables() == 1:
    x = f.univariate_polynomial().roots(multiplicities=False)
    if type(x) == list:
      x = x[0]
    return [{str(f.variables()[0]): x}]
  elif f.nvariables() == 2:
    variables = [str(f.variables()[0]), str(f.variables()[1])]
    sols = []
    if "x1" in variables:
      sols.append(sol[0])
    if "x2" in variables:
      sols.append(sol[1])
    if "x3" in variables:
      sols.append(sol[2])
    if "x4" in variables:
      sols.append(sol[3])

    return coron(f, bounds[0], bounds[1], sol=sols, debug=debug)
  elif f.nvariables() == 3:
    variables = [str(f.variables()[0]), str(f.variables()[1]), str(f.variables()[2])]
    sols = [x for x in sol]
    if "x1" not in variables:
      sols.pop(0)
    elif "x2" not in variables:
      sols.pop(1)
    elif "x3" not in variables:
      sols.pop(2)
    elif "x4" not in variables:
      sols.pop(3)
    
    if method == "coron":
      return coron_trivariate(f, bounds[0], bounds[1], bounds[2], l=0, sol=sols, debug=debug)
    elif method == "bauer-joux":
      return bj07_trivariate(f, bounds[0], bounds[1], bounds[2], opts={"debug":False, "returnFirst": True, "solution": sols, "S": 2})
    else:
      sys.exit("Method `%s` is not implemented. You can use either `coron` or `bauer-joux`." % method)

    return coron(f, bounds[0], bounds[1], sol=sols, debug=debug)
  elif f.nvariables() == 3:
    variables = [str(f.variables()[0]), str(f.variables()[1]), str(f.variables()[2])]
    sols = [x for x in sol]
    if "x1" not in variables:
      sols.pop(0)
    elif "x2" not in variables:
      sols.pop(1)
    elif "x3" not in variables:
      sols.pop(2)
    elif "x4" not in variables:
      sols.pop(3)
    
    if method == "coron":
      return coron_trivariate(f, bounds[0], bounds[1], bounds[2], l=0, sol=sols, debug=debug)
    elif method == "bauer-joux":
      return bj07_trivariate(f, bounds[0], bounds[1], bounds[2], opts={"debug":False, "returnFirst": True, "solution": sols, "S": 2})
    else:
      sys.exit("Method `%s` is not implemented. You can use either `coron` or `bauer-joux`." % method)
    # return bd_trivariate(f, bounds[0], bounds[1], bounds[2], l=2, debug=True)
  else:
    return coron_quadvariate(f, bounds[0], bounds[1], bounds[2], bounds[3], l=0, debug=debug)
