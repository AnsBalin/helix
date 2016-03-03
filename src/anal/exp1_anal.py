
import re 
import os
import math
import glob
import matplotlib.pyplot as plt

def fR(x):

	return re.findall('(R.*)', x)

def fF(x):

	return re.findall('(F.*)', x)

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
				#flt_line = map(float, re.findall(r'-*\d+\.*\d*', line))
				flt_line = [float(i) for i in re.findall(r'-*\d+\.*\d*', line)]
				ply_num = int(flt_line[0])
				data[ply_num][t].append([flt_line[1],flt_line[2],flt_line[3]])





	f.close()
	return data

def anal_fun( Fname, Rname ):

	fdat = readXYZ( Fname )
	rdat = readXYZ( Rname )

	v_sum = lambda v1, v2 : [(v1[0]+v2[0]), (v1[1]+v2[1]), (v1[2]+v2[2])]

	f_tot = []
	r_tot = []
	
	# Force in helix
	#for f_n in fdat[0]:
	#	f_tot.append(reduce( v_sum, f_n ))

	# Positions in polymer
	#for r_n in rdat[1]:
	#	r_tot.append(reduce( v_sum, r_n ))

	#Incomprehensible list comprehension
	hist = [0 for ii in range(0,100)]
	#hists = [ [len( [x for x in [ math.sqrt(y[0]*y[0]+y[1]+y[1]) for y in r_n ] if int(3*x) == i] ) for i in range(0,100)] for r_n in rdat[1]]
	for r_n in rdat[1]:
	#	radii = [ math.sqrt(x[0]*x[0]+x[1]+x[1]) for x in r_n ]
		#radii = map( lambda x: math.sqrt(x[0]*x[0]+x[1]+x[1]), r_n )

		hist = [ (hist[i] + len( [x for x in [ math.sqrt(y[0]*y[0]+y[1]*y[1]) for y in r_n ] if int(3*x) == i] )) for i in range(0,100) ]

	return hist


data_path = "../../exp1/"

files = os.listdir(data_path)
filesR = map(fR, files)
filesF = map(fF, files)

filesR = filter(len, filesR)
filesF = filter(len, filesF)


filesR.sort()
filesF.sort()

narr = [5, 11, 23, 60]
hlar = [5, 15, 25]
xarr = [5, 15, 20]

f_tot = []
sim = 0

hist_n = [[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)]] 

hist_h1 = [[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)]]
hist_x1 = [[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)]]

hist_h2 = [[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)]]
hist_x2 = [[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)],[0.0 for ii in range(0,100)]]

hist_t1 = [0.0 for ii in range(0,100)]
hist_t2 = [0.0 for ii in range(0,100)]


r = [(1.0/3.0)*i for i in range(0,100)]
fig_num=0
for indn, n in enumerate(narr):
	for indh, hl in enumerate(hlar):
		for indx, x in enumerate(xarr):
			
			hist_i  = [0.0 for ii in range(0,100)]
			
			for i in range(0,10):
				print sim
				hist = anal_fun( data_path+filesF[sim][0], data_path+filesR[sim][0] )
				
				hist_i        = [(hist_i[ii] + hist[ii]) for ii in range(0,100)]

				hist_t1       = [(hist_t1[ii]       + (1.0/float(n))*hist[ii]) for ii in range(0,100)]
				hist_t2       = [(hist_t2[ii]                      + hist[ii]) for ii in range(0,100)]
				
				hist_n[indn]  = [(hist_n[indn][ii]                       + hist[ii]) for ii in range(0,100)]

				hist_h1[indh] = [(hist_h1[indh][ii] + (1.0/float(n))*hist[ii]) for ii in range(0,100)]
				hist_h2[indh] = [(hist_h2[indh][ii]                + hist[ii]) for ii in range(0,100)]
				hist_x1[indx] = [(hist_x1[indx][ii] + (1.0/float(n))*hist[ii]) for ii in range(0,100)]
				hist_x2[indx] = [(hist_x2[indx][ii]                + hist[ii]) for ii in range(0,100)]
				
				sim += 1

					#plt.plot( map( lambda x: x[2], f_tot) )
					#plt.plot( [x[2] for x in f_tot] )
			plt.figure(fig_num)
			plt.plot( r[1:], [hist_i[i]/(r[i]*sum(hist_i)) for i in range(1,100)] )
			plt.legend(['n='+str(n)+', l='+str(hl)+', x0='+str(x)])
			plt.savefig('figs/hist'+str(sim)+'.png')
			plt.close(fig_num)
			fig_num += 1
			#plt.show()
			#t = range(0,len(f_tot[0]))



plt.figure(fig_num)
hist_t1 = [hist_t1[i]/(r[i]*sum(hist_t1)) for i in range(1,100)]
plt.plot( r[1:], hist_t1 )
plt.legend(['total'])
plt.savefig('figs/hist_total.png')
plt.close(fig_num)
fig_num += 1


plt.figure(fig_num)
hist_t2 = [hist_t2[i]/(r[i]*sum(hist_t2)) for i in range(1,100)]
plt.plot( r[1:], hist_t2 )
plt.legend(['total2'])
plt.savefig('figs/hist_total2.png')
plt.close(fig_num)
fig_num += 1


plt.figure(fig_num)
hist_n0 = [hist_n[0][i]/(r[i]*sum(hist_n[0])) for i in range(1,100)]
plt.plot( r[1:], hist_n0 )
plt.legend(['n=5'])
plt.savefig('figs/hist_n5.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_n1 = [hist_n[1][i]/(r[i]*sum(hist_n[1])) for i in range(1,100)]
plt.plot( r[1:], hist_n1 )
plt.legend(['n=11'])
plt.savefig('figs/hist_n11.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_n2 = [hist_n[2][i]/(r[i]*sum(hist_n[2])) for i in range(1,100)]
plt.plot( r[1:], hist_n2 )
plt.legend(['n=23'])
plt.savefig('figs/hist_n23.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_n3 = [hist_n[3][i]/(r[i]*sum(hist_n[3])) for i in range(1,100)]
plt.plot( r[1:], hist_n3 )
plt.legend(['n=60'])
plt.savefig('figs/hist_n60.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_h10 = [hist_h1[0][i]/(r[i]*sum(hist_h1[0])) for i in range(1,100)]
plt.plot( r[1:], hist_h10 )
plt.legend(['l=5'])
plt.savefig('figs/hist_l1_5.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_h11 = [hist_h1[1][i]/(r[i]*sum(hist_h1[1])) for i in range(1,100)]
plt.plot( r[1:], hist_h11 )
plt.legend(['l=15'])
plt.savefig('figs/hist_l1_15.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_h12 = [hist_h1[2][i]/(r[i]*sum(hist_h1[2])) for i in range(1,100)]
plt.plot( r[1:], hist_h12 )
plt.legend(['l=25'])
plt.savefig('figs/hist_l1_25.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_h20 = [hist_h2[0][i]/(r[i]*sum(hist_h2[0])) for i in range(1,100)]
plt.plot( r[1:], hist_h20 )
plt.legend(['l=5'])
plt.savefig('figs/hist_l2_5.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_h21 = [hist_h2[1][i]/(r[i]*sum(hist_h2[1])) for i in range(1,100)]
plt.plot( r[1:], hist_h21 )
plt.legend(['l=15'])
plt.savefig('figs/hist_l2_15.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_h22 = [hist_h2[2][i]/(r[i]*sum(hist_h2[2])) for i in range(1,100)]
plt.plot( r[1:], hist_h22 )
plt.legend(['l=25'])
plt.savefig('figs/hist_l2_25.png')
plt.close(fig_num)
fig_num += 1


plt.figure(fig_num)
hist_x10 = [hist_x1[0][i]/(r[i]*sum(hist_x1[0])) for i in range(1,100)]
plt.plot( r[1:], hist_x10 )
plt.legend(['x=5'])
plt.savefig('figs/hist_x1_5.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_x11 = [hist_x1[1][i]/(r[i]*sum(hist_x1[1])) for i in range(1,100)]
plt.plot( r[1:], hist_x11 )
plt.legend(['x=15'])
plt.savefig('figs/hist_x1_15.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_x12 = [hist_x1[2][i]/(r[i]*sum(hist_x1[2])) for i in range(1,100)]
plt.plot( r[1:], hist_x12 )
plt.legend(['x=20'])
plt.savefig('figs/hist_x1_20.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_x20 = [hist_x2[0][i]/(r[i]*sum(hist_x2[0])) for i in range(1,100)]
plt.plot( r[1:], hist_x20 )
plt.legend(['x=5'])
plt.savefig('figs/hist_x2_5.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_x21 = [hist_x2[1][i]/(r[i]*sum(hist_x2[1])) for i in range(1,100)]
plt.plot( r[1:], hist_x21 )
plt.legend(['x=15'])
plt.savefig('figs/hist_x2_15.png')
plt.close(fig_num)
fig_num += 1

plt.figure(fig_num)
hist_x22 = [hist_x2[2][i]/(r[i]*sum(hist_x2[2])) for i in range(1,100)]
plt.plot( r[1:], hist_x22 )
plt.legend(['x=20'])
plt.savefig('figs/hist_x2_20.png')
plt.close(fig_num)
fig_num += 1
