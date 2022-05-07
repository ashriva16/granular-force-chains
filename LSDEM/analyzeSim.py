import numpy as np
import sys
sys.path.insert(1, './preProcess/')
from grain import Grain
import matplotlib.pyplot as plt
from numpy import linalg as LA

def GetPackingFraction(Rve,plotName,title): 
	#get packing fraction 
	grains = Rve.grains
	numParticles = Rve.nGrains
	posRot             = Rve.positions
	posRotEnd          = posRot[Rve.nSteps - 1]
	finalPositions = posRotEnd[:,:2]


	centroid = np.mean(finalPositions,axis=0)
	Ds = []
	PFs = []

	for D in range(100,600,5):

		totalParticleArea = 0

		for particleNum in range(numParticles):
			pos = finalPositions[particleNum]
			grain = grains[particleNum]
			area = grains[particleNum]._area
			if ( (LA.norm(pos - centroid) < D) ):
				totalParticleArea += area         	

		totalArea = np.pi*(D**2)
		PF = totalParticleArea/totalArea
		Ds += [D]
		PFs += [PF]

	plt.plot(Ds,PFs)
	plt.xlabel('Radius from centroid')
	plt.ylabel('Packing density')
	plt.title(title)
	plt.ylim(0.7, 1.0)
	plt.savefig(plotName)
	plt.close()

	return (Ds[25],PFs[25])

	#plt.show()

#get contact number distribution

def getContactNumber(morphDirName,cFileNum,posFileName,plotName,title):
	posFile            = morphDirName + posFileName 
	posRot             = np.loadtxt(posFileName)
	positions          = posRot[:,:2]

	cInfoFile = morphDirName + "c_info/cinfo_" + str(cFileNum) + ".dat"
	cInfo = np.loadtxt(cInfoFile)
	cCount = np.zeros(len(positions)) #ith entry is number of contacts of ith particle 

	centroid = np.mean(positions,axis=0)
	Ds = np.arange(100,400,10)
	pairs = []#list of pairs, to avoid double counting

	for contact in cInfo:  
		p1,p2 = int(contact[0]),int(contact[1]) #particles involved in this contact
		if (not ({p1,p2} in pairs) ): #check if pair already counted 
			position1,position2 = positions[p1],positions[p2]
			centDist1,centDist2 = LA.norm(position1 - centroid),LA.norm(position2 - centroid)

			if (min(centDist1,centDist2) < 200):
				cCount[p1] += 1
				cCount[p2] += 1 
			pairs += [{p1,p2}]


		

	plt.hist(cCount,bins = 100,range=(1,10))
	plt.title(title)
	plt.savefig(plotName)
	plt.close()

	CNav = cCount[np.nonzero(cCount)].mean()
	return CNav

