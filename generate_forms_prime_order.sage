from ast import literal_eval

### GENERATE FORMS WITH MAGMA
magma.load('magma/generate_prime_deg_isogeny.m')

# define the polynomial ring
P.<x1,x2,x3,x4> = PolynomialRing(ZZ)

log_name = "prime_forms.log"
res_name = "prime_forms_res.log"
log = open(log_name, "w")
res_f = open(res_name, "w")

# how many forms per prime
n = 3
print("Generating {} forms per prime...".format(n))

for i in range(1,10+1):
    # gen prime
    p = 0
    while not (p % 4 == 3):
        p = random_prime(2^(100*i), False, 2^(100*i-1))
        log.write("Generated {} bit prime {}\n".format(100*i, p))
        sys.stdout.write("Generated {} bit prime {}\n".format(i, p))
        sys.stdout.flush()
    # gen prime degree
    for j in range(60,120):
        d = random_prime(2^(i * j), False, 2^(i * j - 1))
        log.write("Generated {} bit prime degree {}\n".format(j*i, d))
        sys.stdout.write("Generated {} bit prime degree {}\n".format(j*i, d))
        sys.stdout.flush()
        # call magma to get n forms
        res = magma.eval("RunMeNTimes({},{},{});".format(p, d, n))
        res = res.splitlines()
        forms = []
        sols = []
        # log.write("{}\n".format(str(res)))
        # sys.stdout.write("{}\n".format(str(res)))
        # sys.stdout.flush()
        # continue
        for r in res:
            if r not in "[0]" or len(r) >= 4:
                r = r.replace(',', '')
                # check if sol
                if r[0] == '[':
                    temp = r.replace('  ', ',').replace(' ', ',')
                    lstr = ""
                    while lstr != temp:
                        lstr = temp
                        temp = lstr.replace(',,', ',').replace('[,', '[')
                    sol = literal_eval(lstr)
                    sols.append(sol)
                # otherwise form
                else:
                    form = P(r)
                    forms.append(form)
        log.write("New forms:\n{}\nSols:\n{}\n".format(str(forms),str(sols)))
        sys.stdout.write("New forms:\n{}\nSols:\n{}\n".format(str(forms),str(sols)))
        sys.stdout.flush()
        for data in zip(forms, sols):
            res_f.write("({}, {}, {}, {})\n".format(str(p), str(d), str(data[0]), str(data[1])))

log.close()
res_f.close()
