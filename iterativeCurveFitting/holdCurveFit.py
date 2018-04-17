import numpy as np
from numpy import random
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as pl
import utilities as ut

p0=np.array([1.5,200.,2.,700.,3.,2400.])
hold=[0,1,0,1,0,0]
holdloc=[i for i, x in enumerate(hold) if x]
p1=np.delete(p0,holdloc)
timepoints = np.arange(-200.,4000.,20.)
callFitName= ut.genCallFitName(p0,hold)

y=ut.multiexp(timepoints,p0)+np.random.normal(0, .1, size=len(timepoints))

popt, pcov = curve_fit(ut.multiexp2, timepoints, y, p0=p0)
print 'unheld parameters:', popt

popt, pcov=eval(callFitName)
for i in range(0, len(hold)):
	if hold[i]:
		popt=np.insert(popt,i,p0[i])
print 'held parameters:', popt

yfit=ut.multiexp(timepoints,popt)

pl.plot(timepoints,y,timepoints,yfit)
pl.show()
