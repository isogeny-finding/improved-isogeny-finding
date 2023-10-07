def s(D):
    res = 0
    for i in range(0, D+1):
        res += binomial(i+2, 2)
    return res

def sx(D):
    res = 0
    for i in range(0, D+1+1):
        res += (D + 2 - i) * (i + 1)
    for i in range(0, D+1):
        res += (D + 1 - i) * (i + 1)
    return res

def eps(D):
    return (s(D) / (3 * sx(D) - 2 * s(D))).numerical_approx()

for i in range(1, 100+1):
    sys.stdout.write("({},{})".format(i, round(eps(i),4)))
    sys.stdout.flush()
