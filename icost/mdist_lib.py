def get_gff_intron_probs(gff):

	with open(gff) as fp:
		introns = {}
		total_score = 0
		for line in fp.readlines():
			if line.startswith('#'): continue
			line = line.rstrip()
			line = line.split('\t')
			if len(line) < 8: continue
			if line[2] == 'intron': #and line[6] == '+': why was this here?
				beg = int(line[3])
				end = int(line[4])
				score = line[5]
				intron = (beg, end)
				if score == '.': introns[intron] = 0
				else: 
					if intron not in introns: introns[intron] = 0
					introns[intron] += float(score)
					total_score += float(score)
		
		for i in introns:
			introns[i] = introns[i]/total_score

		return introns

def get_mdist(introns1, introns2):	

	for i in introns1:
		if i not in introns2:
			introns2[i] = 0
	for i in introns2:
		if i not in introns1:
			introns1[i] = 0

	introns1 = dict(sorted(introns1.items(), 
					key=lambda item: item[1], reverse=True))

	dd = 0
	isonum = 0
	for i in introns1:
		beg = i[0]
		end = i[1]
		iprob1 = '{0:.6f}'.format(introns1[i])
		iprob2 = '{0:.6f}'.format(introns2[i])
		#print(f'{introns1[i]}\t{introns2[i]}')
		isonum += 1
		d = introns1[i] - introns2[i]
		dd += abs(d)
		
	return float('{0:.6f}'.format(dd)), isonum

def icost_groups(icost_gffs, wb_gffs):

	mdist_groups = {}
	for ID in icost_gffs:
		for path in icost_gffs[ID]:
			gff1_apc = path 
			gff2_wb = wb_gffs[ID]
			icost = gff1_apc.split('_')[2]
			print('gene ID:', ID)
			print('tested icost:', icost)
			print('apc file path:', gff1_apc)
			print('wb file path:', gff2_wb)
			cap = subprocess.run(f'python3 {program} {gff1_apc} {gff2_wb}', 
				shell=True, capture_output=True)
			print('mdist:', find_mdist(cap))
			print('#####')
			info = [
				{
					'ID': ID, 
					'mdist': find_mdist(cap), 
					'apc_file': gff1_apc.split('/')[-1], 
					'wb_file': gff2_wb.split('/')[-1]
				}
			]
			if icost not in mdist_groups:
				mdist_groups[icost] = info
			else:
				mdist_groups[icost] += info
			#break
		#break

	jsonString = json.dumps(sorted(mdist_groups.items()), indent=4)
	jsonFile = open(args.outdir+args.outfile_name+'.json', 'w')
	jsonFile.write(jsonString)
	jsonFile.close()

