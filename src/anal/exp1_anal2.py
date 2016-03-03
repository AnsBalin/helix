
import re 
import os
import math
import glob
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA
import numpy as np

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

def com( r ):
	return [sum([ x[0]/len(r) for x in r ]), sum([ x[1]/len(r) for x in r ]), sum([ x[2]/len(r) for x in r ])] 

def Rg( r, r_cm ):
	return math.sqrt(sum( [ ((x[0]-r_cm[0])*(x[0]-r_cm[0]) + (x[1]-r_cm[1])*(x[1]-r_cm[1]) + (x[2]-r_cm[2])*(x[2]-r_cm[2])) for x in r ] )/len(r) )

def Re( r ):
	return math.sqrt((r[len(r)-1][0]-r[0][0])*(r[len(r)-1][0]-r[0][0]) + (r[len(r)-1][1]-r[0][1])*(r[len(r)-1][1]-r[0][1]) + (r[len(r)-1][2]-r[0][2])*(r[len(r)-1][2]-r[0][2]) )

def R_PCA( r ):

	myPCA = PCA( np.array(r) )

	return myPCA.fracs[0]/myPCA.fracs[1]

def anal_fun( Fname, Rname ):

	fdat = readXYZ( Fname )
	rdat = readXYZ( Rname )
	radius = []
	
	r_cm_arr = [ com(r) for r in rdat[1] ]
	#radius = [ math.sqrt(r_cm[0]*r_cm[0] + r_cm[1]*r_cm[1]) for r_cm in r_cm_arr ]
	z = [ r_cm[2] for r_cm in r_cm_arr ]
	#R_g = [ Rg(r, r_cm) for (r,r_cm) in zip(rdat[1],r_cm_arr) ] 
	#R_e = [Re(r) for r in rdat[1] ]

	pca_frac = [ R_PCA(r) for r in rdat[1] ]
	return [ z, pca_frac ]


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
fig_num = 0
for indn, n in enumerate(narr):
	f1, axes1 = plt.subplots(3,3,sharex=True,sharey=True)
	#f2, axes2 = plt.subplots(3,3)
	for indh, hl in enumerate(hlar):
		for indx, x in enumerate(xarr):
			
			hist_i  = [0.0 for ii in range(0,100)]
			#plt.figure(fig_num)
			for i in range(0,10):
				print sim
				[z, frac] = anal_fun( data_path+filesF[sim][0], data_path+filesR[sim][0] )
				axes1[indx][indh].plot(z[0::5], frac[0::5])
				
				sim += 1
			#plt.ylim( 0, 40 )
			#plt.xlim( 0, 25 )
			#plt.legend(['n='+str(n)+', l='+str(hl)+', x0='+str(x)])
			f1.savefig('figs/meta/z_frac_'+str(fig_num)+'.png')
			plt.close(f1)
			fig_num += 1
