import re 

def readXYZ( fname ):
	f = open(fname, 'r')

	NUMBER_LINE 	= 0
	COMMENT_LINE	= 1
	DATA_LINE		= 2
	N = 0
	t = -1 					#time
	data = [ [],[] ]
	status = NUMBER_LINE
	for line in f:

		if status == NUMBER_LINE:
			t += 1
			N = int(re.search(r'\d+',line).group())
			status += 1
			data[0].append([])
			data[1].append([])

		elif status == COMMENT_LINE:

			status += 1

		elif status == DATA_LINE:

			N -= 1
			if N == 0:
				status = NUMBER_LINE
			else:
				flt_line = map(float, re.findall(r'-*\d+\.*\d*', line))
				ply_num = int(flt_line[0])
				data[ply_num][t].append([flt_line[1],flt_line[2],flt_line[3]])





	f.close()
	return data


#data = readXYZ("R_211_130.xyz")
#for molecule in data:
#	for time in molecule:
#		for atom in time:
#			print atom 
#		print "\n"
#	print "\n"
