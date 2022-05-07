# -*- coding: utf-8 -*-
"""
Create 2D morph file from surface points
@author: konstantinos
"""
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath
from utilities import smoothHeaviside,sussman

def createMorphs(Rve):
    morphDirName = Rve.morphDir 
    nMorphs      = Rve.nShapes    
 
    # Directories
    ptsDir = morphDirName
    morphDir = morphDirName

    # Grain properties 
    rho = 2.65 #g/pixel^3
    kn = 30000.
    ks = 27000.
    mu = 0.4

    # Misc 
    pad = 2
    eps = 0.5
    plot = False 
    initTimesteps = 5

    for i in range(nMorphs):
        # Read surface points
        ptsFile = ptsDir + str(i) + '.dat'
        pts = np.loadtxt(ptsFile)
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
        
        # Generate the morphology file
        f = open(morphDir + 'grainproperty' + str(i) + '.dat','w')
        # Compute mass (density=1, thickness=1)
        m = rho * smoothHeaviside(-lset, eps).ravel().sum()
        # Compute center of mass
        cx = rho / m * np.dot(smoothHeaviside(-lset, eps), x).sum()
        cy = rho / m * np.dot(smoothHeaviside(-lset, eps).T, y).sum()
        # Compute moment of inertia
        x_dev = np.power(x-cx,2)
        y_dev = np.power(y-cy,2)
        I = rho * np.dot(smoothHeaviside(-lset, eps), x_dev).sum() + \
            rho * np.dot(smoothHeaviside(-lset, eps).T, y_dev).sum()
        # Compute bbox radius
        radius = np.linalg.norm(ptsCM, axis = 1).max()
        pts = ptsCM + np.array([cx,cy])
        # Write to file in the order below
        f.write('%5.3f'% m + '\n')                           # mass
        f.write('%5.3f'% I + '\n')                           # moment of inertia
        f.write('%5.3f' % cx + ' ' + '%5.3f'% cy + '\n')     # center of mass
        f.write('%d' % nPts + '\n')                          # no of boundary points
        f.write(" ".join('%5.6f' % x for x in ptsCM.ravel().tolist()) + '\n') 
        f.write('%5.3f'% radius + '\n')                      # bounding box radius
        f.write('%5.3f' % nX + ' ' + '%5.3f'% nY + '\n')     # lset arr dimensions
        f.write(" ".join(str(x) for x in lset.ravel().tolist()) + '\n') 
        f.write('%5.3f'% Rve.particleProps["kn"] + '\n')                           # kn
        f.write('%5.3f'% Rve.particleProps["ks"] + '\n')                           # ks
        f.write('%5.3f'% Rve.particleProps["mu"] + '\n')                                  # mu
        f.write('%5.3f'% Rve.particleProps["cresN"] + '\n')                        #cresN
        f.write('%5.3f'% Rve.particleProps["cresS"] + '\n')                        #cresS
        f.close()
        # Plot
        if plot:
            fig1 = plt.figure(i)
            ax1 = fig1.add_subplot()
            plt.contour(x, y, lset, 15, linewidths = 0.5, colors = 'k')
            plt.pcolormesh(x, y, lset, cmap = plt.get_cmap('rainbow'))
            plt.colorbar()
            plt.scatter(pts[:,0], pts[:,1],c='r')
            plt.show(block=True)
 
