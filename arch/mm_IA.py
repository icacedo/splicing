s = 'ATGGCTCGATGCTAAGTATGGCTCGATGCT'
k = 3

kmers = {}
mm = {}
for n in range(len(s)-2):
	ctx = s[n:n+2]
	nt = s[n+k-1]
	print(ctx,nt)
	if s[n:n+k] not in kmers:
		kmers[s[n:n+k]]=1
	else:
		kmers[s[n:n+k]]+=1
	if ctx not in mm:
		mm[ctx] = {}
	if nt not in mm[ctx]:	
		mm[ctx][nt] = 0
	mm[ctx][nt]+=1
print(mm)
print(kmers)

