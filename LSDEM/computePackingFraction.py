# -*- coding: utf-8 -*-
"""
Created on Wed June 8 14:55:05 2018

@author: konstantinos
"""

import numpy as np
import scipy.spatial as sp
import matplotlib.pyplot as plt

def getPackingFraction(Rve,plotName,title): 
	
	voronoiVolumes,bdGrains = Rve.getVoronoiVolumes()
	nonBdGrains = list(set(range(Rve.nGrains)).difference(set(bdGrains))) #get grains not on boundary  
	grainVolumes = np.array([grain._area for grain in Rve.grains]) #list of grain areas 
	phi = grainVolumes[nonBdGrains]/voronoiVolumes[nonBdGrains] #packing fraction (grain area/voronoi area) 

	# Compute packing fraction in center RVE of increasing radius
	pos = Rve.positions[-1][nonBdGrains,:2]
	tree = sp.cKDTree(pos) 
	nR = 20
	gRadius = Rve.grains[0]._bBoxRadius
	R = np.linspace(5*gRadius,25*gRadius,nR) #radii averaging 
	phiR = np.zeros(nR)

	for i,r in enumerate(R):
	    idx = tree.query_ball_point(np.mean(pos,axis=0),r) #all grains within a ball of radius r 
	    phiR[i] = np.mean(phi[idx]) #average packing fraction in each ball 

	def movingAverage(a, n=3):
	    ret = np.cumsum(a, dtype=float)
	    ret[n:] = ret[n:] - ret[:-n]
	    return ret[n-1:]/n

	# Converged packing fraction as the one closest to the moving average    
	phiR_ma = movingAverage(phiR,n=3)
	phiR_ma = np.insert(phiR_ma,0,[phiR[0],phiR[0]])
	idxConverged = np.mean(phiR_ma[4:8])#np.argmin(np.abs(phiR-phiR_ma))

	plt.plot(R,phiR,'-bo')
	plt.plot(R,phiR_ma,'-rx')
	plt.ylim(0.7, 0.95)
	plt.title(title)
	plt.savefig(plotName)
	plt.close()
	return (0,idxConverged)
