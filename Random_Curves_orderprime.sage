import sys
import random
rank_firstprime = int(sys.argv[1])
nb_primes= int(sys.argv[2])
pas= int(sys.argv[3])
LL=[]
P=Primes()
for i in range(rank_firstprime,(nb_primes*pas)+rank_firstprime+1,pas):
	p=P.unrank(i)
	while True:
		a = random.randrange(p)
		b = random.randrange(p)
		if ((4* pow(a,3,p))+(27*pow(b,2,p))) % p != 0 :
			E = EllipticCurve(GF(p), [a, b])
			N = E.cardinality()
			if is_prime(N):
				print(p,a,b,N,end=";")
				break

