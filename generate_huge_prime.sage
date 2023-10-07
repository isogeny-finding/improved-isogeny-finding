from ast import literal_eval

load('coron_trivariate.sage')
load('bj07_trivariate.sage')

### GENERATE FORMS WITH MAGMA
magma.load('magma/random_forms.m')

# define the polynomial ring
P.<x1,x2,x3,x4> = PolynomialRing(ZZ)

# log file
log_name = "huge_primes.log"
# how many bits max, 100*k bits for 1..50
k = 50
# how many primes
n_primes = 1
# how many forms per prime
n_forms = 10

log = open(log_name, "w")

for i in range(1, k+1):
    log.write("Prime bit length: {}, will generate {} primes.\n".format(100*i, n_primes))
    for j in range(1, n_primes+1):
        # generate random prime
        p = 0
        while not (p % 4 == 3):
            log.write("Generating prime...\n")
            p = random_prime(2^(100*i), False, 2^(100*i-2))
        p_len = len(bin(p))-2
        print("\n=== Prime: ({}) {}\n".format(p_len, p))
        log.write("=== Prime: ({}) {}\n".format(p_len, p))
        # steps to take in form, was set to prime bit length
        # steps = p.bit_length()
        steps = p_len
        # call magma to get n forms
        res = magma.eval("RunMeNTimes({},{},{});".format(p, steps, n_forms))
        res = res.splitlines()
        Q = []
        for r in res:
            if r not in "[0]" or len(r) > 4:
                r = r.replace(',', '')
                form = P(r)
                Q.append(form)
        log.write("[\n")
        for ii in range(n_forms):
            print("Form: {}".format(str(Q[ii])))
            log.write("{},\n".format(str(Q[ii])))
        log.write("]\n")

print("Done")
log.close()
