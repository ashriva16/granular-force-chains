#Use genetic algorithms to produce a clone in 2D
#to work at multiple length scales
#To run, import this file (GE.py) and call makeParticles(inputDir)
#inputDir should contain input file ("inputFile.dat"), which is of form:
#aspectRatio \n Mu_roundness \n Std_roundness \n Mu_circularity \n Std_circularity \n numParticles \n

import numpy as np 
import math
import random
import scipy
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy import interpolate 
from scipy.interpolate import splprep, splev
import copy
import pylab

from deap import tools
from deap import base 
from deap import creator 

import matplotlib.path as mpltPath
from utilities import smoothHeaviside,sussman
import numpy.linalg as la
from scipy.spatial import Delaunay
import scipy.stats as stats


#Parameters chosen by user
Std_r = 0.5 #standard deviation of radius mutation 
Std_th = 0.1 #Standard deviation of theta mutation 
circularityLB = 0.8 #lower bound on circularity 
radBound     = 10 #Only record curvature at points with ROC below this 
Plot         = False #Plot grains

init_stepsize = 1 #initial step size
popSize = 20 #size of population 
MUTP = 0.1 #probability of mutation 
NGEN = 0 #Number of generations
CXPB = 0.2 #cross-over probability
weight_r = 0 #weight on match radius, 0 < weight_r < 1
weight_c = 1 - weight_r
Iterations = 2
Subdivisions = 2

numOutputPoints = 100 #How many points to output

def cartToPol(xl,yl): #convert xl (x list) and yl (y list) to rl and Thetal 
	length = len(xl)
	rl,Thetal = np.zeros(length),np.zeros(length)
	for i in range(length):
		xv,yv = xl[i],yl[i]
		rl[i] = (xv**2 + yv**2)**0.5
		Thv = math.atan2(float(yv),xv)
		if (Thv < 0): Thv = 2*math.pi + Thv
		Thetal[i] = Thv
	return rl,Thetal

#convert polar coordinates (r[i],Theta[i]) to (x,y)
def polToCart(r,Theta,i):
	(rval,theta) = r[i],Theta[i]
	x = rval*math.cos(theta)
	y = rval*math.sin(theta)
	return x,y

#Make element list and flag list (list which states which vertices correspond to which length scales), and theta list, up to It (current iteration)
def makeLists(It):
	if (It):
		points_above = 8*(Subdivisions**(It - 1))
	else: points_above = 0

	points = 8*(Subdivisions**It)
	e_length = points
	element_list = np.ones((e_length,4))*-1
	e_pos = 0 #position in element list
	for i in range(points):
		weight = 1
		element_list[e_pos] = [i, (i -  1) % points, (i + 1) % points, weight]
		e_pos += 1

	return element_list


#calculate curvature given radii list and index into radii list
def calc_curvature(x,y,last,next,i): 
	(x0,y0) = x[last],y[last]
	(x1,y1) = x[i],y[i]
	(x2,y2) = x[next],y[next]
	x = [x0,x1,x2]
	y = [y0,y1,y2]
	tck,uout = splprep([x,y],s=0.,per=False,k=2)
	xd1,yd1 = splev(uout,tck,der = 1) #x'(t),y'(t)
	xd2,yd2 = splev(uout,tck,der = 2) #x''(t),y''(t)
	x1_1,y1_1 = xd1[1],yd1[1] #x1'(p1),y1'(p1)
	x2_1,y2_1 = xd2[1],yd2[1] #x2'(p1),y2'(p1)
	curvature = (x1_1*y2_1 - y1_1*x2_1)/( (x1_1**2 + y1_1**2)**(3./2) )
	return curvature

#calculate spline length given radii list and index into radii list
def calc_perimeter(x,y): 
	nPoints = len(x)

	tck,u = splprep([x,y],s=0.,per=False,k=3)
	u_new = np.linspace(u.min(), u.max(), 100)
	xd1,yd1 = splev(u_new,tck,der = 1) #x'(t),y'(t)
	arcVals = np.sqrt(xd1**2 + yd1**2) #integrate along arc length to get perimeter
	du      = u_new[1] - u_new[0]
	perimeter = scipy.integrate.simps(y = arcVals, x = u_new, dx = du)

	return perimeter

#mutation function for particle

def mutatePart(particle,element_list):
	points = len(particle[0])
	for element in element_list:
		point = int(element[0])
		if (random.random() < MUTP): 
			rval = particle[0][point]
			Thetaval = particle[1][point]
			particle[0][point] = np.random.normal(rval,Std_r)
			particle[1][point] = np.random.normal(Thetaval, Std_th)
	return particle

#Calculate the first and last index of the element list for this level
def calc_num_e(It):
	if not (It):
		return 0,8
	else:
		e_start = 8*(Subdivisions**(It - 1))
		e_end = 8*(Subdivisions**(It))
		return e_start, e_end

#Get lset from particle 
def getLsetPolar(particle):
	x,y = getCartVals(particle)
	tck, u = splprep([x,y], u=None, s=0.0, per=1) 
	u_new = np.linspace(u.min(), u.max(), 50)
	x_new, y_new = splev(u_new, tck, der=0) #First, interpolate more points
	lset = getLset(x_new,y_new)
	return lset


def getLset(X,Y): #get the level set for a set of points (from Kostas)
	pad = 2
	initTimesteps = 5


	pts = np.column_stack((X,Y))
	ptsCM = pts - np.mean(pts,axis=0)
	nPts = len(pts)
    # Create grid to evaluate level set
	xMin, yMin = np.min(ptsCM,axis=0)
	xMax, yMax = np.max(ptsCM,axis=0)
	cm = [-xMin+pad,-yMin+pad]
	pts = ptsCM + cm
	nX = int(np.ceil(xMax-xMin+2*pad))
	nY = int(np.ceil(yMax-yMin+2*pad))
	x = np.arange(0,nX)
	y = np.arange(0,nY)
	xx,yy = np.meshgrid(x,y)
    # Evaluate signed distance on the grid
	path = mpltPath.Path(pts)
	lset = np.zeros(xx.shape)
	for j in range(nY):
		for k in range(nX):
			xy = [xx[j,k],yy[j,k]]
			dist = la.norm(xy-pts,axis=1)
			idx = np.argmin(dist)
			lset[j,k] = dist[idx]
			inside = path.contains_points([xy])
			lset[j,k] *= -1 if inside else 1
    # Reinitialize level set
	for j in range(initTimesteps):
		lset = sussman(lset,0.1)

	return lset


def getLsetProps(x,y,lset): #get the maximum insribed sphere in a particle & area from lset
	maxRad = abs(np.amin(lset))
	area   = np.count_nonzero(lset <= 0)
	return maxRad, area

#evaluate particle fitness
def evalPart(particle, element_list, It, generation):
	cost = 0
	length = len(element_list[:,0])
	particleNum = particle.particleNumber
	radiusCurvatures = np.zeros(length) #list of radii of curvatures 

	x,y = getCartVals(particle)
	tck, u = splprep([x,y], u=None, s=0.0, per=1) 
	u_new = np.linspace(u.min(), u.max(), numOutputPoints)
	x, y = splev(u_new, tck, der=0)#points generated from smoothe spline 
	
	tck, u = splprep([x,y], u=None, s=0.0, per=1) #re-fit spline to smooth points 
	u_new = np.linspace(u.min(), u.max(), length + 1)
	x, y = splev(u_new, tck, der=0) #Get new points 
	

	for e_num in range(length): #get curvature costs
		element = element_list[e_num]
		i = int(element[0])
		last = int(element[1])
		next = int(element[2])
		weight = int(element[3])  #element weight
		cval = calc_curvature(x,y, last, next, i) #current curvature at this point
		radiusCurvatures[e_num] = abs(1./cval)

	lset = particle.lset
	maxRadius,area = getLsetProps(x,y,lset) #get maximum radius of an inscribed circle & particle area

	radiusCurvaturesNormed = abs(radiusCurvatures/maxRadius) #Radius of curvatures normed by max radius
	radiusCurvaturesNormedFiltered = np.array([radiusCurvaturesNormed[i]  if radiusCurvaturesNormed[i] < radBound else radBound for i in range(length)]) #radius of curvatures less than 1 - major surface features
	
	roundness = np.mean(radiusCurvaturesNormedFiltered)
	perimeter = calc_perimeter(x,y) 
	circularity = 2*np.sqrt(np.pi*area)/perimeter
	particle.roundness = roundness 	
	particle.RCs = radiusCurvaturesNormedFiltered

	if (It == 0): 
		cost = ( (circularity - circularityTarget)/circularity)**2 
	else: 
		cost =  ((roundness - roundnessTarget)/roundness)**2
		cost += len(radiusCurvaturesNormedFiltered[radiusCurvaturesNormedFiltered < 0.4])*100 #penalize large curvatures
 

	print (generation,circularity,circularityTarget,roundness,roundnessTarget, particleNum,cost)
	if (cost < 0.001): print (radiusCurvaturesNormedFiltered)
	#print particle.lset
	return (cost,)

def getCartVals(part): #Get list of cartesian values in order given from r,theta
	length = len(part[0])
	r, Theta = part[0], part[1]
	x = np.zeros(length + 1)
	y = np.zeros(length + 1)
	for i in range(length):
		xv,yv = polToCart(r,Theta,i)
		x[i],y[i] = xv,yv
	x[length], y[length] = x[0],y[0]
	return x,y

def plot(part, It,count): #plot a 2d particle from list of r values)
	plt.figure()
	x,y = getCartVals(part)
	tck, u = splprep([x,y], u=None, s=0.0, per=1, k = 3) 
	u_new = np.linspace(u.min(), u.max(), 10000)
	x_new, y_new = splev(u_new, tck, der=0)
	plt.plot(x_new, y_new, label = "1st Iteration")

	for i in range(len(x) - 1):
		plt.text(x[i], y[i], str(part[2][i]), fontsize=12)


	plt.axes().set_aspect('equal')
	plt.show()
	plt.close()

#output points into file
def outputPoints(part,particleNum, morphDir):
	x,y = getCartVals(part)
	tck, u = splprep([x,y], u=None, s=0.0, per=1) 
	u_new = np.linspace(u.min(), u.max(), numOutputPoints)
	x, y = splev(u_new, tck, der=0)#points generated from spine
	pts = np.column_stack((x,y))
	outName = morphDir + str(particleNum) + ".dat"
	np.savetxt(fname = outName, X = pts, fmt="%f")

	tck, u = splprep([x,y], u=None, s=0.0, per=1) #re-fit spline
	u_new = np.linspace(u.min(), u.max(), 17) #get new points
	x, y = splev(u_new, tck, der=0) #Get new points on smoothed spline

#re-initialize population after each iteration
def re_init_pop(pop,It,newPart):
	points = 8*(Subdivisions ** (It + 1)) #number of points in new iteration
	for partnum in range(popSize): #iterate through population
		part = copy.deepcopy(newPart)
		x,y = getCartVals(part)
		tck, u = splprep([x,y], u=None, s=0.0, per = 1) #fit spline to current particle
		u_new = np.linspace(u.min(), u.max(), points, endpoint = False) 
		x_new, y_new = splev(u_new, tck, der = 0) 
		r, Theta = cartToPol(x_new, y_new)
		pop[partnum][0] = r
		pop[partnum][1] = Theta

	return pop

#check that ind doesn't have any non-physical (high curvature) features at a lower length scale 
def isPhysical(ind,element_list): 
	evalPart(ind, element_list, 2, 3)
	RCs = ind.RCs[1:-1]
	print (RCs)
	return (np.amin(RCs) > 0.1)

#initialize genetic algorithm
def make_GE():
	global toolbox
	#initialize DEAP genetic algorithm population
	points = 8
	creator.create("FitnessMin", base.Fitness, weights = (-1.0,))
	creator.create("Particle", list, fitness = creator.FitnessMin, stepsize = init_stepsize)
	toolbox = base.Toolbox()
	toolbox.register("rand_float", random.randint,1,1)
	toolbox.register("pointList", tools.initRepeat, list, toolbox.rand_float, points)
	toolbox.register("particle", tools.initRepeat, creator.Particle, toolbox.pointList, 2)
	toolbox.register("population", tools.initRepeat, list, toolbox.particle, popSize)
	toolbox.register("evaluate",evalPart)
	toolbox.register("mate",tools.cxTwoPoint)
	toolbox.register("select", tools.selTournament, tournsize = 2)

#main program - calls all others
def clone(particleNum, aspectRatio, Mu_roundness, Std_roundness, Mu_circularity, Std_circularity, Area, morphDir):
	global roundnessTarget
	global circularityTarget

	make_GE()
	pop = toolbox.population() #initialize population
	Min_prin = np.sqrt(Area/(np.pi*aspectRatio))
	Max_prin = Min_prin*aspectRatio
	
	for part in pop: #Some extra initialization not handled by DEAP framework
		Theta = np.linspace(0,2*math.pi,8, endpoint = False)
		r = (Min_prin*Max_prin)/np.sqrt( (Min_prin*np.cos(Theta))**2 + (Max_prin*np.sin(Theta))**2 )
		part[0] = r
		part[1] = Theta
		part.particleNumber = particleNum
		part.lset = getLsetPolar(part)

	count = 0
	element_list = makeLists(0)
	for It in range(Iterations):
		

		roundnessTarget = np.random.normal(Mu_roundness,Std_roundness)  #get random target curvature
		circularityTarget = stats.truncnorm.rvs(
    		(circularityLB - Mu_circularity) / Std_circularity, (1 - Mu_circularity) / Std_circularity, loc=Mu_circularity, scale=Std_circularity, size = 1)[0]
		
		fitnesses = [toolbox.evaluate(pop[i],element_list,It,0) for i in range(len(pop)) ] #get fitnesses of each individual in population
		for ind,fit in zip(pop,fitnesses): #assign fitness values to individuals in population
			ind.fitness.values = fit

		for g in range(NGEN): #begin evolutionary process through generations
			print("-- Generation %i --" % g)
			#Select next generation of individuals through tournament 
			offspring = toolbox.select(pop,popSize)

			#Clone the selected individuals
			offspring = list(map(toolbox.clone, offspring))

			#Apply crossover and mutation on the offspring 

			for child1, child2 in zip(offspring[::2], offspring[1::2]):
				if random.random() < CXPB: 
					child1, child2 = toolbox.mate(child1, child2)
					del child1.fitness.values 
					del child2.fitness.values
					child1.lset = getLsetPolar(ind) #get new lset for mutant
					child2.lset = getLsetPolar(ind) #get new lset for mutant



			#Now randomly mutate

			for mutant in offspring:
				if (random.random() < MUTP):
					mutant = mutatePart(mutant,element_list)
					del mutant.fitness.values
					mutant.lset = getLsetPolar(ind) #get new lset for mutant

			#calculate the fitness of any offsprings w/ invalid fitness
			invalid_inds = [ind for ind in offspring if not ind.fitness.valid]
			fitnesses = [toolbox.evaluate(invalid_ind,element_list,It,g) for invalid_ind in invalid_inds]
			for ind,fit in zip(invalid_inds, fitnesses):
				ind.fitness.values = fit


			#Set the population to be equal to the new offspring, then repeat!
			pop[:]  = offspring
			fitnesses = [toolbox.evaluate(pop[i],element_list,It,g) for i in range(len(pop)) ] 
			#if (g%25 == 0):
			#	plot(pop[0], It, count)
				#print min(fitnesses)	
			count += 1
			
			if ( (np.amin(fitnesses) < 10**-3 and It == 1) or (np.amin(fitnesses) < 10**-6 and It == 0) ): break  
			if (g == (NGEN - 1)): 
				return 0	

		if (Plot): plot(pop[0], It, count) #set up next length scale 
		bestInd = pop[np.argmin(fitnesses)]
		pop = re_init_pop(pop,It,bestInd)
		element_list = makeLists(It + 1)

	if (not isPhysical(pop[0],element_list)): return 0 
	outputPoints(pop[0], particleNum, morphDir)
	count += 1
	return 1 

#program initialize parameters and genetic algorithm
def makeParticles(Rve,particleNum):
	Area = math.pi*(10)**2
	print (Rve.morphDir)

	(aspectRatio, Mu_roundness, Std_roundness, Mu_circularity, Std_circularity, numParticles, morphDir) = \
		Rve.aspectRatio,Rve.Mu_roundness,Rve.Std_roundness,\
		Rve.Mu_circularity,Rve.Std_circularity,int(Rve.nShapes), Rve.morphDir

	while (not clone(particleNum, aspectRatio, Mu_roundness, Std_roundness, Mu_circularity, Std_circularity, Area, morphDir)):
		pass


