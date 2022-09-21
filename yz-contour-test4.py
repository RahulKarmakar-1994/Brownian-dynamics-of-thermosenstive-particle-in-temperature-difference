import numpy as np
import math
import os
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
cwd = os.getcwd()
sigma = 1                  #sigma for water. Use as length scale.
binwidth = 0.1
Lx = 39.5                         #box sides in sigma = 1 units
L = 18.0

nbins = int(math.ceil(L/binwidth +2 ))

bin_vol = binwidth*binwidth
nf = 101

f = open(cwd +'/../'+'Steadystate-Coordinates1.dat' , 'r')

xlw = [-9.875]#Layer - CL1
xlh = [-12.5]

#time = [1600000,1700000,1800000,1900000,2000000]
time = [1800000]
#time = [2300000]
for t in time:
	tr = t + 50000+1
	nm = str(int(t*0.0001))
	t_stop = np.arange(t,tr,500)
	#print(t_stop)
	k = len(t_stop) - 1
	mm = np.zeros((len(xlw),nbins,nbins))
	count = 0
	while True :
		s = f.readline().split()
		if not s: break 
		if int(s[4]) in t_stop :
			x = float(s[0]) - Lx*round(float(s[0])/Lx)
			y = float(s[1]) - L*round(float(s[1])/L)
			z = float(s[2]) - L*round(float(s[2])/L)
			d = int(s[3])
			for layer in range(len(xlw)):
				if x <= xlw[layer] and x >= xlh[layer] :
					binx = int(math.floor((y)/binwidth))
					biny = int(math.floor((z)/binwidth))
					mm[layer][binx,biny] = mm[layer][binx,biny] + 1.0
					count = count + 1
		if int(s[4]) > t_stop[k] : break		
	#print(count)			
	mm = mm/(bin_vol*nf)
#	heat = max(max(i) for i in mm[layer])
#	mm = mm/heat
	#print(mm)
	xx = np.zeros(nbins)
	yy = np.zeros(nbins)
	for i in range(nbins):
		xx[i] = (i+0.5)*binwidth*sigma
		yy[i] = (i+0.5)*binwidth*sigma
	X, Y = np.meshgrid(xx, yy)
	levels = np.arange(0,20,0.2)
	for layer in range(len(xlw)):
		plt.rcParams["font.weight"] = "bold"
		#plt.rcParams["axes.labelweight"] = "bold"
		fig, ax = plt.subplots(figsize=(8,6), dpi=100)
		ax.spines["top"].set_linewidth(5)
		ax.spines["left"].set_linewidth(5)
		ax.spines["right"].set_linewidth(5)
		ax.spines["bottom"].set_linewidth(5)
		#cs = ax.contourf(X, Y, mm[layer],levels, cmap='gist_gray')
		cs = ax.contourf(X, Y, mm[layer],levels, cmap='gist_heat_r')
		#levels = np.arange(0,0.3,0.001)
		#plt.contourf(X, Y, d2,levels, cmap='Set1')
		bounds = [0,4,8,12,16,20]
		cbar = fig.colorbar(cs,orientation = 'vertical' , pad = 0.25,ticks=bounds)
		cbar.ax.tick_params(labelsize=30)
		cbar.ax.tick_params(labelsize=30)
		ax.set_xticks(np.arange(0, 18, 9))
		ax.set_yticks(np.arange(0, 18, 9))
		ax.xaxis.set_major_locator(MultipleLocator(9))
		ax.yaxis.set_major_locator(MultipleLocator(9))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

# For the minor ticks, use no labels; default NullFormatter.
		ax.xaxis.set_minor_locator(MultipleLocator(4.5))
		ax.yaxis.set_minor_locator(MultipleLocator(4.5))
		ax.tick_params(axis="both", direction="in", length=15, width=4, color="black")
		ax.tick_params(axis="x", labelsize=30, labelrotation=0)
		ax.tick_params(axis="y", labelsize=30, labelrotation=0)
		ax.tick_params(axis = 'x',which='minor',  direction="in",length=10, width=4,color='black')
		ax.tick_params(axis = 'y',which='minor',  direction="in",length=10, width=4,color='black')
		ax.set_xlabel('Z' , size = 35 ,fontweight='bold')
		ax.xaxis.set_label_coords(0.5, -0.11)
		ax.set_ylabel("Y", size = 35,fontweight='bold')
		ax.yaxis.set_label_coords(-0.07, 0.5)
		fig.savefig('Contourmap-test2-colour3-'+nm+'2.png', bbox_inches='tight',dpi=600)
		fig.savefig('Contourmap-test2-colour3-'+nm+'-pres2.png', bbox_inches='tight',dpi=100)
		#plt.show()
		#plt.clf()
	
