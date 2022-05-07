# -*- coding: utf-8 -*-
"""

@author: konstantinos
@date: November, 2016
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from grain import Grain
import helper
from utilities import *

def plotSim(Rve,figDir,vidName,lX,lY,makePics = 1,video = 1,forces = 0):
    print ("forces",forces)
    #get objects
    morphDir = Rve.morphDir
    morphFile = Rve.morphFile
    nGrains            = Rve.nGrains
    posRot             = Rve.positions
    posRot0            = Rve.startPositions
    nSteps             = Rve.nSteps
    morphID = Rve.morph 
    print ("figDir",figDir)

    rho = 2.65 #g/pixel^3   

    if not os.path.exists(figDir):
        os.mkdir(figDir)

    # Read grain morphology data

    # Instantiate grains
    grains = Rve.grains

    if not (makePics):
        nSteps = 0

    print ("nSteps",nSteps)
    for step in range(nSteps-30,nSteps,1):
        # Set up the figure
        fig, ax = plt.subplots(figsize = (5,5))
        ax.set_xlim(0,lX)
        ax.set_ylim(0,lY)
        ax.autoscale_view()
                
        # Update grain positions
        posRotStep = posRot[step]
        for n in range(nGrains):
            grains[n].updateXY(posRotStep[n])
        
        # Collect in patches
        patches = []
        for n in range(nGrains):
            poly = Polygon(grains[n]._newPoints, True) 
            patches.append(poly)

        # Setup up patchCollection 
        pCol = PatchCollection(patches, facecolors='dimgray',edgecolors='black', lw=0.1)
        ax.add_collection(pCol)

        if forces and step == nSteps - 1:
            cInfoFiles = ["cinfo_0.dat"]
            cInfoFinal = cInfoFiles[-1] #get filename of final cInfo 
            print (Rve.cInfoDir + cInfoFinal) 
            cinfo =  np.loadtxt(Rve.cInfoDir + cInfoFinal)

            # Visualize contacts 
            cLen = 30.; cWid = 2    
            forceChains = computeForceChains(cinfo,cLen,cWid,'b')
            ax.add_collection(forceChains)  

        #plt.axis('off')
        plt.savefig(figDir + '/step_' + str(step) + '.png', format='png', dpi =  500 )
        #plt.show(block=True)
        plt.close(fig)

    # Create Video
    # if (video): 
    #     string = "ffmpeg -start_number 0 -i %sstep_" % figDir + "%d.png -y -vcodec mpeg4 " + vidName
    #     print (string)
    #     os.system(string)
    
    return grains
