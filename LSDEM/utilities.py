import numpy as np

def smoothHeaviside(X, eps):
    cond1 = (np.absolute(X) <= eps).astype(int)
    cond2 = (X > eps).astype(int)
    res = cond1*0.5*(1 + X/eps + np.sin(np.pi*X/eps)/np.pi) + cond2
    return res

def sussman(D, dt):
    # forward/backward differences
    a = D - np.roll(D, 1, axis=1)
    b = np.roll(D, -1, axis=1) - D
    c = D - np.roll(D, -1, axis=0)
    d = np.roll(D, 1, axis=0) - D

    a_p = np.clip(a, 0, np.inf)
    a_n = np.clip(a, -np.inf, 0)
    b_p = np.clip(b, 0, np.inf)
    b_n = np.clip(b, -np.inf, 0)
    c_p = np.clip(c, 0, np.inf)
    c_n = np.clip(c, -np.inf, 0)
    d_p = np.clip(d, 0, np.inf)
    d_n = np.clip(d, -np.inf, 0)

    a_p[a < 0] = 0
    a_n[a > 0] = 0
    b_p[b < 0] = 0
    b_n[b > 0] = 0
    c_p[c < 0] = 0
    c_n[c > 0] = 0
    d_p[d < 0] = 0
    d_n[d > 0] = 0

    dD = np.zeros_like(D)
    D_neg_ind = np.flatnonzero(D < 0)
    D_pos_ind = np.flatnonzero(D > 0)

    dD.flat[D_pos_ind] = np.sqrt(
        np.max(np.concatenate(
            ([a_p.flat[D_pos_ind]**2], [b_n.flat[D_pos_ind]**2])), axis=0) +
        np.max(np.concatenate(
            ([c_p.flat[D_pos_ind]**2], [d_n.flat[D_pos_ind]**2])), axis=0)) - 1
    dD.flat[D_neg_ind] = np.sqrt(
        np.max(np.concatenate(
            ([a_n.flat[D_neg_ind]**2], [b_p.flat[D_neg_ind]**2])), axis=0) +
        np.max(np.concatenate(
            ([c_n.flat[D_neg_ind]**2], [d_p.flat[D_neg_ind]**2])), axis=0)) - 1

    D = D - dt * sussman_sign(D) * dD
    return D

def sussman_sign(D):
    return D / np.sqrt(D**2 + 1)

import numpy as np
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def computeCumulativeForceChains(cinfo, cLen=20, cWid=0.5, color='black'):
    
    if cinfo.ndim > 1:
        uniqVal, uniqInv = np.unique(cinfo[:,0:2],return_inverse=True,axis=0)
        cCum = np.zeros((len(uniqVal),9))
        # Compute cumulative contacts
        for i,uv in enumerate(uniqVal):
            dupIdx = np.where(uniqInv == i)[0]
            cCum[i,0:2] = cinfo[dupIdx[0],0:2]
            cCum[i,2:4] = np.sum(cinfo[dupIdx,2:4], axis=0)
            cCum[i,4:6] = np.mean(cinfo[dupIdx,4:6], axis=0)
            cCum[i,6:8] = np.mean(cinfo[dupIdx,6:8], axis=0)
            cCum[i,8] = len(dupIdx)
        if cCum.shape[0] > 1:
            cLocs = cCum[:,6:8]
            cForces = cCum[:,2:4]
            cForceMags = np.linalg.norm(cForces,axis=1)
            # Get forces higher than avg and normalize them    
            cForceMax = np.percentile(cForceMags,95)
            mask = cForceMags > cForceMax
            nMask = len(np.where(mask)[0])
            cForces[mask] = [cForces[mask][i]/cForceMags[mask][i]*cForceMax for i in range(nMask)]
            cForcesNormalized = cForces/cForceMax        
            cForcesNormalizedMags = np.linalg.norm(cForcesNormalized,axis=1)
            
            # Create collection of lines representing force chains
            lines = []
            for i in range(len(cLocs)):
                [x1, y1] = cLocs[i] - cForcesNormalized[i]*cLen/2
                [x2, y2] = cLocs[i] + cForcesNormalized[i]*cLen/2
                lines.append([(x1,y1),(x2,y2)])  
            
            # Collect all to a LineCollection
            lwd = cForcesNormalizedMags * cWid 
            lCol = LineCollection(lines, color=color, linewidths=lwd)

        else:
            cLocs = cCum[0][6:8]
            cForces = cCum[0][2:4]
            cForces = cForces/np.linalg.norm(cForces)
            [x1, y1] = cLocs - cForces*cLen/2
            [x2, y2] = cLocs + cForces*cLen/2
            lines = [[(x1,y1),(x2,y2)]]  
            lwd = cWid 
            lCol = LineCollection(lines, color=color, linewidths=lwd)
            
    else:
        cLocs = cinfo[6:8]
        cForces = cinfo[2:4]
        cForces = cForces/np.linalg.norm(cForces)
        [x1, y1] = cLocs - cForces*cLen/2
        [x2, y2] = cLocs + cForces*cLen/2
        lines = [[(x1,y1),(x2,y2)]]  
        lwd = cWid 
        lCol = LineCollection(lines, color=color, linewidths=lwd)

    return lCol
    
def computeForceChains(cinfo, cLen=20, cWid=0.5, color='black'):
    
    if cinfo.ndim > 1:
        cLocs = cinfo[:,6:8]
        cForces = cinfo[:,2:4]
        cForceMags = np.linalg.norm(cForces,axis=1)
        # Get forces higher than avg and normalize them    
        cForceMax = np.percentile(cForceMags,95)
        mask = cForceMags > cForceMax #bool array of super high forces 
        nMask = len(np.where(mask)[0]) #number of grains satisfying condition above 
        print (mask,nMask,cForces[mask])
        cForces[mask] = [cForces[mask][i]/cForceMags[mask][i]*cForceMax for i in range(nMask)] #set biggest forces = 95% 
        cForcesNormalized = cForces/cForceMax #normalize all forces by 95% force 
        cForcesNormalizedMags = np.linalg.norm(cForcesNormalized,axis=1) #get magnitude of normalized forces
        
        # Create collection of lines representing force chains
        lines = []
        for i in range(len(cLocs)):
            [x1, y1] = cLocs[i] - cForcesNormalized[i]*cLen/2
            [x2, y2] = cLocs[i] + cForcesNormalized[i]*cLen/2
            lines.append([(x1,y1),(x2,y2)])  
        
        # Collect all to a LineCollection
        lwd = cForcesNormalizedMags * cWid 
        lCol = LineCollection(lines, color=color, linewidths=lwd)
    else:
        lCol = LineCollection([], color=color, linewidths=cWid)

    return lCol

