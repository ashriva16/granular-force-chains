# -*- coding: utf-8 -*-
"""
Grain class 

@author: Konstantinos Karapiperis
"""
import numpy as np
class Grain():
    
    def __init__(self, morphFile):
        self._nPoints = 0 
        self._points = []
        self.readMorphFile(morphFile)
    
    def readMorphFile(self, morphFile):
        '''
        Reads the grain property file and stores variables
        '''
        f = open(morphFile, 'r')
        lines = f.readlines()
        self._area = float(lines[0])/2.65
        self._nPoints = int(lines[3])
        self._bBoxRadius = float(lines[5])
        self._points = np.array([float(val) for val in lines[4].split()])
        self._points = np.reshape(self._points,(self._nPoints,2))
        self._newPoints  = np.copy(self._points)
        self._shrinkDistance = 0.5*self.getPerimeter()/len(self._points)     
        self.cmLset = np.array([float(val) for val in lines[2].split()])
        self.nLset = np.array([int(float(val)) for val in lines[6].split()])
        self._xdim, self._ydim = self.nLset
        lset = np.array([float(val) for val in lines[7].split()])                               
        self.lset = np.reshape(lset,self.nLset,order='F').ravel(order='C')        
        f.close()

    def shrinkGrid(self, shrinkDistance = 'None'):
        '''
        Shrinks surface point grid towards the interior by shrinkDistance
        '''
        if shrinkDistance == 'None':
            shrinkDistance = self._shrinkDistance
        
        for i in range(self._nPoints):
            normal = self.getNormal(self._points[i])
            self._points[i] += shrinkDistance * normal

    def getNormal(self,point):
        '''
        Finds normal to level set at a given point
        '''
        x,y = point
        xf = int(np.floor(x))
        yf = int(np.floor(y))
        xc = int(np.ceil(x))
        yc = int(np.ceil(y))
        b1 = self.getGridValue(xf,yf) 
        b2 = self.getGridValue(xc,yf)-b1
        b3 = self.getGridValue(xf,yc)-b1
        b4 = -b2 - self.getGridValue(xf,yc) + self.getGridValue(xc,yc)  
        dx = x-xf
        dy = y-yf
        gradient = [b2+b4*dy, b3+b4*dx]
        gradient /= np.linalg.norm(gradient)  
        return gradient

    def getGridValue(self,x,y):
        '''
        Finds value of level set at a given point
        '''
        return self.lset[y*self._xdim + x]

    def getUpdatedPoints(self, posRot):
        '''
        Updates grain configuration by rotating and moving the _points
        '''
        # Translate to the origin, rotate and translate back
        theta = posRot[2]
        cm = posRot[0:2]
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta), np.cos(theta)]])
        points =  self._points.dot(R.T) + cm
        return points


    def updateXY(self,posRotNew):
        self._newPoints = self.getUpdatedPoints(posRotNew)


    def getPerimeter(self):
        '''
        Computes perimeter
        '''
        perimeter = 0
        for i in range(len(self._points)):
            pt1 = self._points[i]
            pt2 = self._points[(i+1)%len(self._points)]
            perimeter += np.linalg.norm(pt2-pt1)
        return perimeter