# -*- coding: utf-8 -*-
"""
@author: Konstantinos Karapiperis
and Robert Buarque de Macedo 
"""
import GE #my libraries 
import create_morphology 
import plotConfiguration
import computePackingFraction
import analyzeSim
import helper
from grain import Grain
from copy import deepcopy


import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial as sp
import networkx as nx
from numpy import linalg as LA
import numpy.ma as ma 
import os
import pathlib
import pickle
from joblib import Parallel, delayed
import multiprocessing

class rve():

    #default values for rve instances
    def __init__(self,morphDirName,mainDirName):  
        self.aspectRatio = 1.0 #AR = max axis/min axis of starting ellipse
        self.Mu_roundness = 1.0  #Average of roundness targets
        self.Std_roundness = 0.01 #standard deviation of roundness targets
        self.Mu_circularity = 1.0 #same, but for circularity
        self.Std_circularity = 0.01 
        self.trial = 0
        self.nShapes  = 1 #number of morph geometries to be created
        self.nGrains = 900 #number of grains in simulation 
        self.morphDir = morphDirName + "/"*(morphDirName[-1] != "/") #Add "/" on end if it is not there 
        self.mainDir = mainDirName + "/"*(mainDirName[-1] != "/")

        self.cInfoDir = self.mainDir + "c_info/" 
        self.morphFile  =  "morphIDs.dat"  #folder which maps morphID (index) 
        self.posFile = "positionsOut.dat" #Output files for simulation 
        self.velFile = "velocitiesOut.dat"
        self.startPosFile = "positions.dat" #input files for simulation 
        self.startVelFile = "velocities.dat"
        self.paramFile = "param_file.txt"
        self.shearHistFileIn = "shearIn.dat" #shear history input file name, set as "None" if none
        self.shearHistFileOut = "shearHist.dat"#shear hist output file name 

        self.numParticlesRow = 15 #pluviation parameters for putting particles in a grid
        self.numParticlesCol = self.nGrains/self.numParticlesRow 
        self.grainRadiusVert = 40 #distances between particles in grid
        self.grainRadiusHorz = 100./3
        self.startSpot = [2*self.grainRadiusHorz,2*self.grainRadiusVert] #bottom left of grid

        self.nTot = 500000 #simulation parameters: number of steps
        self.nOut = 10000 #how often data is output
        self.stopStep = 1900000 #step to stop shaking
        self.A = 0.0 #Amplitude of shake
        self.T = 0.0 #time period of shake 
        self.topWallPos = (4 + self.numParticlesCol)*self.grainRadiusVert #position of top wall 
        self.rightWallPos = 600 #position of right wall (left and bottom are at x=0/y=0)
        self.Force = 0.0  #Force on top wall 
        self.dropDist = 50. #for pluviation mode 
        self.simDir = "LSDEM2D/" #directory where LSDEM simulation files are 
        self.scriptRunFile = self.simDir + "scriptrun.sh" #filename of executable for running program
        self.randomMult = 200 #random number which is maximum starting velocity in a direction 

        self.particleProps = {} #material parameters of the particles   
        self.particleProps["mu"] = 0.0
        self.particleProps["kn"] = 1e11 
        self.particleProps["ks"] = 1e11 
        self.particleProps["cresN"] = 0.4
        self.particleProps["cresS"] = 0.5

        self.makeDir(self.morphDir)
        self.makeDir(self.mainDir)
        self.makeDir(self.cInfoDir)

#=============================================================================================================
    #load data from a simulation into self
    #After a simulation is run, this function MUST be called in order to 
    #perform any plotting/analysis on the simulation results 

    def getDataFromFiles(self):
        morphDir = self.morphDir
        mainDir = self.mainDir 
        self.morph = np.loadtxt(mainDir + self.morphFile, dtype=int)[1:] #map from particle number to ID
        try:  
            positions = np.loadtxt(mainDir + self.posFile) #get output positions 
            velocities = np.loadtxt(mainDir + self.velFile) #'' velocities #  
        except: #in case a simulation was stopped randomly, skip final time step
            positions = np.genfromtxt(mainDir + self.posFile,skip_footer=1)
            velocities = np.genfromtxt(mainDir + self.velFile,skip_footer=1)
            raise Exception("Unfinished file")

        self.nSteps = int(len(positions)//self.nGrains)
        self.positions = np.array(np.split(positions[:self.nSteps*self.nGrains,:], self.nSteps)) #add dimension for timestep 
        self.velocities = np.array(np.split(velocities[:self.nSteps*self.nGrains,:], self.nSteps))
        self.startPositions = np.loadtxt(mainDir + self.startPosFile)

        self.grains = np.empty(self.nGrains, dtype=object)             
        for n in range(self.nGrains):
            propFile =  morphDir + "grainproperty" + str(int(self.morph[n]))+ ".dat"
            self.grains[n] = Grain(propFile)

        print('Found ', self.nGrains, ' grains')
        print('Found ', self.nSteps, ' steps')
#=============================================================================================================  
    #create new directory dirName if it does not exist 

    def makeDir(self,dirName):
        if not os.path.exists(dirName):
            pathlib.Path("./" + dirName).mkdir(parents=True, exist_ok=True)
#=============================================================================================================
    #Show attribute values of rve instance 

    def showVals(self):
        print( ' '.join("%s: %s \n" % item for item in vars(self).items()) )

#=============================================================================================================
    #make particles using genetic algorithm and convert into level set using create_morphology

    def createParticles(self):
        num_cores = multiprocessing.cpu_count() #use GE to make particles in parallel 
        #Parallel(n_jobs=num_cores)(delayed(GE.makeParticles)(self,i) for i in np.arange(self.nShapes)) 
        for i in np.arange(self.nShapes): GE.makeParticles(self,i)
        create_morphology.createMorphs(self) #convert them into morph files
#=============================================================================================================
    #make initial positions for a grid of particles, given number of particles in row/col,
    #spacing between particles in both directions (grainRadiusVert/grainRadiusHorz), and the origin of the grid 
    #also randomly assign particles to each grid point 

    def makeInitPositions(self):
        numParticlesRow, numParticlesCol,grainRadiusVert,grainRadiusHorz,startSpot = \
            self.numParticlesRow, self.numParticlesCol,self.grainRadiusVert,self.grainRadiusHorz,self.startSpot

        if (self.nGrains != self.numParticlesRow*self.numParticlesCol): #check consistency 
            print ("error: number of particles in row * col != nGrains")
            return

        numParticles = int(numParticlesRow*numParticlesCol)
        positions = np.zeros((numParticles,3))
        count = 0
        for ynum in np.arange(numParticlesCol): #assemble particles in a grid [startSpot_i,startSpot_i + grainRadius*numGrains_i]
            ypos = startSpot[1] + grainRadiusVert*ynum  
            for xnum in np.arange(numParticlesRow):
                xpos = startSpot[0] + grainRadiusHorz*xnum
                positions[count,:2] = [xpos,ypos]
                positions[count,2] = (2*np.pi - -2*np.pi) * np.random.random_sample() + -2*np.pi
                count += 1
        print ("maxYpos",ypos)
        self.startPositions = positions
        np.savetxt(fname = self.mainDir + self.startPosFile, X = positions, fmt = "%f")   
        velocities = np.append(np.zeros((numParticles, 2)), np.zeros([numParticles, 1]), 1)
        # velocities = np.append(400*np.random.random((numParticles,2)) - 200,np.zeros([numParticles,1]),1) #randomly sample starting velocities between -200&200
        np.savetxt(fname = self.mainDir + self.startVelFile, X = velocities, fmt = "%f")

        #make map from grain number to morph 
        Q,R = self.nGrains//self.nShapes,self.nGrains%self.nShapes 
        IDs = np.zeros(self.nGrains)
        for q in np.arange(Q): IDs[q*self.nShapes:(q+1)*self.nShapes] = np.random.permutation(self.nShapes)  
        if (R != 0): IDs[-R:] = np.random.permutation(R)      
        IDs = np.append(self.nGrains, IDs)
        np.savetxt(fname = self.mainDir + self.morphFile, X = IDs, fmt = "%d", delimiter = "\n") #make morphFile

#=============================================================================================================
    #create input parameter file and run sim
    #if run is False, an execution bash script is created as file self.scriptRunFile
    #executing this bash file runs LSDEM with the parameter file containing all simulation information 
    #if run is True, LSDEM is run directly from here. 

    def executeDEM(self,pluviate=False,run=False):
        fname = self.mainDir + self.paramFile
        mainDirName = "../" + self.mainDir
        morphDirName = "../" + self.morphDir
        cInfoDirName = "../" + self.cInfoDir

        parameter_string = """morphIDs = %s
inPos   = %s
inVel   = %s
morphDir = %s
outPos  = %s
outVel  = %s
cInfoDir = %s
shearHistFileIn = %s
shearHistFileOut = %s
nTot    = %d
nOut    = %d
stopStep = %d 
pluviate = %d
numRow   = %d
A       = %f
T       = %f 
Force   = %f
topWallPos = %f
dropDist = %f 
rightWallPos = %f"""%(mainDirName +  self.morphFile,mainDirName + self.startPosFile,mainDirName + self.startVelFile,
            morphDirName, mainDirName + self.posFile, mainDirName + self.velFile, cInfoDirName,mainDirName + self.shearHistFileIn, mainDirName + self.shearHistFileOut, self.nTot, self.nOut,self.stopStep, 
            pluviate, self.numParticlesRow,self.A,self.T,self.Force, self.topWallPos, self.dropDist, self.rightWallPos, 
            )
        text_file = open(fname, "w")

        text_file.write(parameter_string)

        text_file.close()

        runCommand = "mpirun -np 4 ./Main2 " + mainDirName + self.paramFile + ";" #generate script for running
        if (run): os.system("cd " + self.simDir + ";" + runCommand ) #run the simulation
        else: #Don't run, just append run command to a scriptRun file 
            runCommand = "sbatch -A geomechanics HPCMainTest.sh " + mainDirName + self.paramFile
            if os.path.exists(self.scriptRunFile): append_write = 'a' # append if already exists
            else: append_write = 'w' # make a new file if not
            with open(self.scriptRunFile, append_write) as filetowrite:
                filetowrite.write("%s\n" % runCommand)

        return 1
#=============================================================================================================
    #plot domain given X and Y box length (lX/lY), location of where static figures are saved (figDir),
    #name of video to be created from figures (vidName), whether or not figures are made (makePics), and whether 
    #or not a video is made (makeVideos) 
 
    def plot(self,lX,lY,figDir,vidName,makePics = 0,makeVideos = 0):
        positionsEquil = self.positions[-1]#save the last positions of the particles for the next sim 
        equilFileName = "positionsEquil.dat" #save equilibration
        np.savetxt(fname = self.mainDir + equilFileName, X = positionsEquil, fmt = "%f") 
        os.system("rm " + figDir + "*")
        self.makeDir(figDir)
        plotConfiguration.plotSim(self,figDir,vidName,lX,lY,makePics,makeVideos,forces=1)
#=============================================================================================================
    #get packing fraction of a simulation, given directory where packing fraction plot will be (plotDir),
    #name of packing fraction plot (plotName), title in plot (title) and if Kostas' method is used (bool kostas) 

    def getPackingFraction(self,plotDir,plotName,title,kostas = False):
        self.makeDir(plotDir)
        if (kostas): (D,PF) = computePackingFraction.getPackingFraction(self,plotName,title)
        else: (D,PF) = analyzeSim.GetPackingFraction(self,plotName,title)
        return (D,PF)

#=============================================================================================================
    #save rve for future use as 'name' 

    def saveObj(self, name ):
        self.makeDir(name)
        with open(name + '.pkl', 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

#=============================================================================================================
    #get grain volume (area) from morph file 

    def getGrainVolume(self): 
        f = open(self.morphDir+'grainproperty'+str(self.morph[0])+'.dat')
        fr = f.read().splitlines()
        grainVolume = float(fr[0])
        f.close()
        return grainVolume


#=============================================================================================================
    # Compute set Voronoi volumes for all grains at given step 

    def getVoronoiVolumes(self):
        # Initialize
        voronoiVolumes = np.zeros(self.nGrains)
       
        # Grab positions/rotations
        posRot = self.positions[-1]
        pos = posRot[:,:2]

        #get a grain and shrink it
        grain = deepcopy(self.grains[0])
        grain.shrinkGrid()
        
        # For neighbor computations
        tree = sp.cKDTree(pos) 
        bBoxRadius = grain._bBoxRadius
        
        # Keep track of grains with cell boundary intersection
        bdGrains = []
        
        # Save Voronoi cells for visualization ? (useful also for debugging)
        vizCell = 0 #1

  
        # Get updated surface points for all grains
        pts = []
        for n in range(self.nGrains):
            pts.append(grain.getUpdatedPoints(posRot[n]))
        
        # For each grain run tesselation on a subregion centered around the grain 
        for n in range(self.nGrains):
            
            # Find surrounding grains 
            neighbors = tree.query_ball_point(pos[n],4*bBoxRadius)
            
            # Initialize points for Voronoi computation
            cellPts = []
            
            # Add center grain points (plus the cm)
            cellPts.append(pos[n])
            for pt in pts[n]:
                cellPts.append(pt)

            for p in neighbors:
                
                # Center grain already considered
                if p == n: continue
            
                # Restrict to points close to the surface points of the center grain
                keptIdx = []
                for idx in range(len(pts[p])):
                    distFromCenterGrain = np.linalg.norm(pts[p][idx]-pts[n],axis=1)
                    if distFromCenterGrain.min() > 2*bBoxRadius: continue
                    keptIdx.append(idx)
                    
                # Continue if surface points of the neighbor are far from center grain
                if len(keptIdx) == 0: continue
                
                ptsKept = pts[p][keptIdx]
                cellPts.append(pos[p])
                for pt in ptsKept:
                    cellPts.append(pt)
            
            # Run Voronoi computation 
            cellPts = np.array(cellPts)
            vor = sp.Voronoi(cellPts, qhull_options='QJ')  
            subCells = [vor.regions[vor.point_region[i]] for i in range(len(pts[n])+1)]
                            
            # Check that subcells of the center grain don't intersect boundaries
            bdCheck = True
            for i in range(len(pts[n])+1):
                if any(sc == -1 for sc in subCells[i]):
                    bdCheck = False
                    break
            
            if bdCheck == False:
                bdGrains.append(n)
                #print('Boundary grain!', n)
            else:
                # Add volume contributions for subcells of center grain
                vorVertices = vor.vertices
                for sc in subCells[:len(pts[n])+1]:
                    voronoiVolumes[n] += polyArea(vorVertices[sc])

            # For voronoi visualization only
            if vizCell == False: continue
            if bdCheck == False:
                sp.voronoi_plot_2d(vor)
                plt.plot(pts[n][:,0],pts[n][:,1],lw=2,c='r',alpha=0.5)
                plt.show(block=True)
        
        return voronoiVolumes,bdGrains

#=============================================================================================================
    #Compute cumulative contact data from multiple single pair contact data
    #Return: cPair,cForces,cNormals,cLocs

    def getCumulativeContacts(self,cData):
       # Run numpy's unique to the contact pair data
        uniqVal, uniqInv = np.unique(cData[:,0:2],return_inverse=True,axis=0)
        
        # Array holding cumulative values (cPairs, cForces, cNormals, cLocs)
        cDataCum = np.zeros((len(uniqVal),8))

        # Compute cumulative contacts
        for i,uv in enumerate(uniqVal):
            dupIdx = np.where(uniqInv == i)[0]
            cDataCum[i,0:2] = cData[dupIdx[0],0:2]
            cDataCum[i,2:4] = np.sum(cData[dupIdx,2:4], axis=0)
            cDataCum[i,4:6] = np.mean(cData[dupIdx,4:6], axis=0)
            cDataCum[i,6:8] = np.mean(cData[dupIdx,6:8], axis=0)

        return cDataCum

#=============================================================================================================
   #Save coordination number for each grain 

    def getAvgCoordinationNumber(self,title):
        # Read cinfo
        cInfoFiles = helper.getSortedFileList(self.cInfoDir)
        cInfoFinal = "cinfo_0.dat"#cInfoFiles[-1] #get filename of final cInfo  
        cAll =  np.loadtxt(self.cInfoDir + cInfoFinal)
        cAll = self.getCumulativeContacts(cAll)
        cPairs = cAll[:,:2].astype(int)
        centroid = cAll[:,6:].mean(axis = 0)
        endPositions = self.positions[-1]
        # Use networkx binary graph
	# Overkill here but it is an opportunity to learn how to use this nice library
        G = nx.Graph()
        # Add edges
        for i,c in enumerate(cPairs):
            G.add_edge(c[0],c[1])

        # Find coordination number
        Z = np.zeros(self.nGrains)
        for i in range(self.nGrains):
            # Continue if rattler
            if i not in G: continue
            if (LA.norm(endPositions[i,:2] - centroid) > 250): continue
            Z[i] = len(G[i])

        ZMasked = ma.masked_array(Z, Z < 1)
        return ZMasked.mean()

#=============================================================================================================
#Shoelace computation of polygon area
 
def polyArea(pts):
    lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
    area = 0.5*abs(sum(x1*y2-x2*y1 for x1,y1,x2,y2 in lines))
    return area
