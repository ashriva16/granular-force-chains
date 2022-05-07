#generate/analyze simple circles in square domain 

import helper
import plotConfiguration
import analyzeSim
import computePackingFraction
from rve import rve
import numpy as np

# PATH = "/home/ankit/Desktop/My_local_work/project_temp/Granular_project/results/LSDEM/"
PATH=""
def run():
    trials = 3
    for trial in np.arange(2,trials):
        # directory where morph files are stored
        morphDirName = PATH+"MLstuff/trial_%d/" % trial
        # directory where all relevant information to specific simulation is stored
        mainDirName = PATH+"MLstuff/trial_%d/" % trial

        Rve = rve(morphDirName,mainDirName) #create instance of rve class for this simulation 
        Rve.nShapes = 2 #number of different unique morphs 
        Rve.nGrains = 121  #number of grains in the experiment 
        Rve.simDir = "LSDEM2D/" #directory where simulation will take place 
        Rve.nTot = 3000 #total timesteps 
        Rve.nOut = 50  #how often data is output from LSDEM 
        Rve.shearHistFileIn = "None" #shear histories input file name. 'None' if there's none 
        Rve.shearHistFileOut = "shearHistories.dat" #shear histories output file name
        Rve.trial = trial

        Rve.numParticlesRow = 11 #pluviation parameters for putting particles in a grid
        Rve.numParticlesCol = Rve.nGrains/Rve.numParticlesRow 
        Rve.grainRadiusVert = 20 #distances between particles in grid (vertical )
        Rve.grainRadiusHorz = 20 #'' (horizontal)
        # Rve.startSpot = [Rve.grainRadiusHorz,Rve.grainRadiusVert] #bottom left of grid
        Rve.startSpot = [10.1,
                         10.1]  # bottom left of grid
        Rve.rightWallPos = 220

        Rve.showVals() #show current attributes 
        Rve.createParticles() #create particles given default roundness/circularity distribution moments
        Rve.makeInitPositions() #set up initial positions/velocities 
        Rve.executeDEM(pluviate=0,run=1) #run simulation. If run == 1, runs locally. If run == 0, creates a scriptrun.sh file that, if executes, submits job on supercomputer  
        Rve.saveObj(PATH+"MLstuff/rvesML/rveML_%d"%trial) #save this rve in directory rveTest

def plot():
    print("-------------------Plotting----------------------------")
    rveDir = PATH+"MLstuff/rvesML/"
    rveFileList = helper.getSortedFileList(rveDir)#get list of saved rve's
    print (rveFileList,rveDir)
    for rveFile in rveFileList:
        
        Rve = helper.loadObj(rveDir + rveFile) #load this rve 
        Rve.getDataFromFiles() #extract data from simulation 
        outputDir = "MLstuff/trial_%d/"%Rve.trial
        Rve.plot(Rve.rightWallPos,Rve.rightWallPos,outputDir + "figs/",outputDir + "vid.mp4",makePics = 1, makeVideos = 1) #plot simulation 
        #except: continue
run()
plot()
