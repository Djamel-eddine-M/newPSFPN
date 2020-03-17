import sys
P = Primes()
interval=int(sys.argv[1])
nb_max=int(sys.argv[2])
for i in range(10,(nb_max*interval)+10,interval):
	print(P[i], end=" ")

