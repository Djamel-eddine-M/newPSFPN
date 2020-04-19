import sys
import random
p = int(sys.argv[1])
while True:
	a = random.randrange(p)
	b = random.randrange(p)
	if ((4* pow(a,3,p))+(27*pow(b,2,p))) % p != 0 :
		E = EllipticCurve(GF(p), [a, b])
		N = E.cardinality()
		if is_prime(N):
			print(a,b,N)
			break

