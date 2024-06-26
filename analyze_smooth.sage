load('find_roots.sage')
load('smooth_forms.sage')

def shortish_iso(method="coron"):
  P.<x1, x2, x3, x4> = PolynomialRing(ZZ)

  examples = get_forms()
 
  S = spline([(58, 3), (6, 2), (101, 2), (121, 2), (141, 6), (151, 14), (161, 17), (166, 19), (171, 21), (181, 26), (201, 37), (211, 42), (251, 50), (331, 40), (351, 49), (391, 70), (396, 72)])

  debug = False
  results = []

  variables_num = 3

  known_variables = [0, 1]
  if variables_num == 3:
    known_variables = [0]
  elif variables_num == 4:
    known_variables = []

  for ind, example in enumerate(examples):
    if ind != 3: # choose which form to run
      continue
    
    p = example["p"]
    forms = example["forms"]

    print("=================== new example ", ind)
    print("p bit length: ", p.bit_length())
 
    for forms1 in forms:
      failureShortSameNum = 0
      failureDependentNum = 0
      successNum = 0
      deg_bit_lens = []

      for ind, f in enumerate(forms1["forms_with_same_ideal_norm"]):
          sol = forms1["sols"][ind]
          iso_deg = forms1["iso_deg"]
          f = f - iso_deg

          bit_len = int(iso_deg).bit_length()
          deg_bit_lens.append(bit_len)

          if bit_len < 150:
            b = 2**13
          else:
            b = 2**ceil(S(bit_len))
          bounds = (b, b, b, b)

          constant_polys = []
          for v in known_variables:
            if v == 0:
              f = f(x1 = sol[0])
              constant_polys.append(x1 - sol[0])
            elif v == 1:
              f = f(x2 = sol[1])
              constant_polys.append(x2 - sol[1])
            elif v == 2:
              f = f(x3 = sol[2])
              constant_polys.append(x3 - sol[2])
            elif v == 3:
              f = f(x4 = sol[3])
              constant_polys.append(x4 - sol[3])

          if debug:
            print("==========")
            print("f: ", f)
            print("known variables: ", known_variables)
            print("solution: ", sol)

          assert f(x1 = sol[0], x2 = sol[1], x3 = sol[2], x4 = sol[3]) == 0

          roots, err = small_roots(f, constant_polys, bounds, p, sol, m=2, method=method, debug=debug)
          if debug:
            print("roots: ", roots)
            print("errMsg", err)

          if len(roots) > 0:
            successNum += 1
          else:
            if err == "short same":
              failureShortSameNum += 1
            else:
              failureDependentNum += 1

          if len(roots) == 0:
            continue

          for root in roots:
            if debug:
              print("---")
              print("root: ", root)
            g = f
            if 0 in known_variables:
              g = f(x1 = sol[0])
            if 1 in known_variables:
              g = f(x2 = sol[1])
            if 2 in known_variables:
              g = f(x3 = sol[2])
            if 3 in known_variables:
              g = f(x4 = sol[3])
            
            root_dict = {}
            for k, v in root.items():
              root_dict[str(k)] = v 
            
            keys = root_dict.keys()

            if "x1" in keys:
              g = g(x1 = root_dict["x1"])
            if "x2" in keys:
              g = g(x2 = root_dict["x2"])
            if "x3" in keys:
              g = g(x3 = root_dict["x3"])
            if "x4" in keys:
              g = g(x4 = root_dict["x4"])

            if debug:
              print("eval: ", g)
            assert g == 0

      print("=========")
      print("success: ", successNum)
      print("failure short same: ", failureShortSameNum)
      print("failure dependent: ", failureDependentNum)
      print(deg_bit_lens)


if __name__ == '__main__':
  if len(sys.argv) > 1:
    method = sys.argv[1]
    shortish_iso(method)
  else:
    shortish_iso()
