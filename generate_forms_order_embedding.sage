load('coron_trivariate.sage')
load('bj07_trivariate.sage')
### GENERATE FORMS WITH MAGMA
magma.load('magma/order_embedding.m')

data = [
    {
        "p": 2^250 * 68 - 1,
        "Q": [],
        "e": 17,
        "E": 21,
        "stats": {},
    },
    # {
    #     # sqi sign prime 3923
    #     "p": 2^65 * 5^2 * 7 * 11 * 19 * 29^2 * 37^2 * 47 * 197 * 263 * 281 * 461 * 521 * 3923 * 62731 * 96362257 * 3924006112952623 - 1,
    #     "Q": [],
    #     "e": 17,
    #     "E": 22,
    #     "stats": {},
    # },
    # {
    #     # sqi sign prime 6983
    #     "p": 2^33 * 5^21 * 7^2 * 11 * 31 * 83 * 107 * 137 * 751 * 827 * 3691 * 4019 * 6983 * 517434778561 * 26602537156291 - 1,
    #     "Q": [],
    #     "e": 17,
    #     "E": 22,
    #     "stats": {},
    # },
    {
        # SIKEp434
        "p": 2^216 * 3^137 - 1,
        "Q": [],
        "e": 25,
        "E": 34,
        "stats": {},
    },
    # {
    #     # SIKEp503
    #     "p": 2^250 * 3^159 - 1,
    #     "Q": [],
    #     "e": 25,
    #     "E": 45,
    #     "stats": {},
    # },
    {
        # SIKEp610
        "p": 2^305 * 3^192 - 1,
        "Q": [],
        "e": 35,
        "E": 45,
        "stats": {},
    },
    # {
    #     # SIKEp751
    #     "p": 2^372 * 3^239 - 1,
    #     "Q": [],
    #     "e": 40,
    #     "E": 80,
    #     "stats": {},
    # },
]

# define the polynomial ring
P.<x1,x2,x3,x4> = PolynomialRing(ZZ)

# how many forms per prime
n = 100
print("Generating {} forms per prime...".format(n))

for d in data:
    p = d["p"]
    print("Prime:", p)
    # steps to take in form, was set to prime bit length
    # steps = p.bit_length()
    steps = len(bin(p))-2
    # call magma to get n forms
    res = magma.eval("RunMeNTimes({},{},{});".format(p, steps, n))
    res = res.splitlines()
    key = str(p)
    for r in res:
        if r not in "[0]":
            r = r.replace(',', '')
            form = P(r)
            d["Q"].append(form)
    assert n == len(d["Q"])

print("Generated forms.")

### EXPERIMENTS

for d in data:
    print("Starting experiments for prime {}:".format(d["p"]))
    for i in range(d["e"], d["E"]+1):
        success = 0
        failure = 0
        for Q in d["Q"]:
            x0 = randint(0, 2^i)
            y0 = randint(0, 2^i)
            z0 = randint(0, 2^i)
            f = Q - Q(x0, y0, z0, 0)
            b = 2^(i+1)
            roots, _ = coron_trivariate(f, b, b, b, l=0, debug=False)
            if roots:
                success += 1
            else:
                failure += 1
        d["stats"][str(i)] = (success, failure)
        print("Exponent:", i, "(S, F):", success, failure, "Success rate:", success / (success + failure))

print("Writing log.")
fd = open("order_embedding.log", "w")
fd.write(str(data))
fd.close()

print("Finished")
