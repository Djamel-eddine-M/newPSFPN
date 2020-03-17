import sys
E=EllipticCurve(GF(sys.argv[1]),[sys.argv[2],sys.argv[3]])
A=E.cardinality()   
print(A)
