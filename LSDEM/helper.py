#helper functions for runDEM that are not in a class
import os 
import pathlib
import numpy as np
import natsort
import pickle5 as pickle


def loadObj(name):
	with open(name, 'rb') as f:
		print (f)
		return pickle.load(f)

def getSortedFileList(Dir):
	files = [f for f in os.listdir(Dir) if os.path.isfile(os.path.join(Dir, f)) and f[-3:] == "pkl" ]
	return natsort.natsorted(files)

