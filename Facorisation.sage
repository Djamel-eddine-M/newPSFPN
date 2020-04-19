import sys
N=int(sys.argv[1])
r=factor(N)
res=list(r)
for i in range(len(res)):
	for j in range(2):
		print(res[i][j], end=" ")

