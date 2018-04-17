
""" 
.. module: FitResFFTPixel0.py
   :platform: Windows 10 Python 2.7
.. moduleauthor:: Scott Ellis <skellis@berkeley.edu> 
A comprehensive bundle of mathematical functions and routines for fitting experimental data. 
Provides methods for
 
    - exponentials
    - generate dummy data simulating time resolved stimulated Raman signals,
    - plotting functions for trouble shooting
    - iterative exponential fitting
    - iterative fourrier transform
..
"""
import numpy as np
from numpy import random
import scipy as sp
from scipy.optimize import curve_fit
import time
import fitfuncs as ff



#A little function for converting wavenumbers to femtoseconds. Also functions in reverse fs to wavenumbers (cm-1).
def wn2fs(wn):
	return 1/(2.99*10**-5*wn)


def trimData(data,timepoints,st):
	desiredLoc=(int(st)-int(timepoints[0]))/(int(timepoints[1])-int(timepoints[0]))
	dataTrim=np.delete(data,np.arange(desiredLoc),axis=1)
	timeTrim=np.delete(timepoints,np.arange(desiredLoc))
	print 'The data has been trimmed at '+str(st)+ ' fs. timepoints has been shortened from '+str(len(timepoints))+' points to '+str(len(timeTrim))+' points.'
	return dataTrim,timeTrim


def generateDummyData(shiftx,timepoints,FSRSCoords,ISRSCoords,expParams,noise):
#FSRSCoords=np.array([amp,freq,fwhm...])
	dummyData=np.zeros((len(shiftx),len(timepoints)),dtype=int)
	t,w=np.meshgrid(timepoints,shiftx)
	
	if len(FSRSCoords)%5!=0:
		print 'wrong number of parameters in FSRSCoords. See comments.'
		return
#FSRSCoords=np.array([amp,phase,freq,fwhm,tauPeak...])
	for i in range(0, len(FSRSCoords)/5):	
		#populationSig=FSRSCoords[0]*FSRSCoords[2]/2/3.1415/((w-FSRSCoords[1])+(FSRSCoords[2]/2)^2)*np.heaviside(t,0.5)*np.exp(-t/FSRSCoords[3])+.0001
		#for now each peak will decay as a single exponential.
		populationSig=FSRSCoords[5*i+0]*(np.cos(FSRSCoords[5*i+1])+np.sin(FSRSCoords[5*i+1])*(w-FSRSCoords[5*i+2]))*FSRSCoords[5*i+3]**2/((w-FSRSCoords[5*i+2])**2+FSRSCoords[5*i+3]**2)*np.heaviside(t,0)*np.exp(-t/FSRSCoords[5*i+4])+.0001
		dummyData=np.add(dummyData,populationSig)
	
	if len(expParams)%2!=0:
		print 'wrong number of parameters in ExpParams. See comments.'
		return
#ExpParams=np.array([amp,tau1...])
	for i in range(0, len(expParams)/2):
		background=expParams[2*i+0]*np.exp(-t/expParams[2*i+1])*np.heaviside(t,0)
		dummyData=np.add(dummyData,background)
	
	if len(ISRSCoords)%3!=0:
		print 'wrong number of parameters in ISRSCoords. See comments.'
		return
#ISRSCoords=np.array([amp,freq,taudephasing])
	for i in range(0, len(ISRSCoords)/3):
		coherenceBackground=ISRSCoords[3*i+0]*np.cos(2*3.14152*t/wn2fs(ISRSCoords[3*i+1]))*np.exp(-t/ISRSCoords[3*i+2])*np.heaviside(t,0)
		dummyData=np.add(dummyData,coherenceBackground)
	
	if 1:
		for i in range(0, len(FSRSCoords)/5):
			for j in range(0, len(ISRSCoords)/3):	
				esFSRSSig=ISRSCoords[3*j+0]*FSRSCoords[5*i+0]/10*(np.cos(FSRSCoords[5*i+1]+t/wn2fs(ISRSCoords[3*j+1]))+np.sin(FSRSCoords[5*i+1]+t/wn2fs(ISRSCoords[3*j+1]))*(w-FSRSCoords[5*i+2]))*FSRSCoords[5*i+3]**2/((w-FSRSCoords[5*i+2])**2+FSRSCoords[5*i+3]**2)*np.heaviside(t,0)*np.exp(-t/FSRSCoords[5*i+4])+.0001
				dummyData=np.add(dummyData,esFSRSSig)
			#
		#
	#

	noiseSig=np.random.normal(0.0,noise,(len(shiftx),len(timepoints)))
	dummyData=np.add(dummyData,noiseSig)
	#plotFSRSContour(populationSig,t,w)
	return dummyData,t,w


def formatBounds(guessCoefs,hold,minAmp,maxAmp,minTau,maxTau):
	holdloc=[i for i, x in enumerate(hold) if x]
	p1=np.delete(guessCoefs,holdloc)
	lowerbound=[0.]*len(guessCoefs)
	upperbound=[0.]*len(guessCoefs)
	lowerbound[0::2]=[minAmp]*int(len(guessCoefs)/2)
	lowerbound[1::2]=[minTau]*int(len(guessCoefs)/2)
	upperbound[0::2]=[maxAmp]*int(len(guessCoefs)/2)
	upperbound[1::2]=[maxTau]*int(len(guessCoefs)/2)
	lowerbound=np.delete(lowerbound,holdloc)
	upperbound=np.delete(upperbound,holdloc)
	bounds=(lowerbound,upperbound)
	return bounds

#note method hacks curve_fit so that it can include the functionality for a hold array along with an arbitrary number of exponentials determined by .
def genCallFitName(p0,hold,bounds=0):
	counter=0
	callFitName='curve_fit(lambda timepoints, *p: ('
	for i in range(0, len(hold)):
		if hold[i] and i%2==0:
			callFitName+=str(p0[i])
		elif hold[i] and i%2==1:
			callFitName+='*np.exp(-timepoints/'+str(p0[i])+')'
		elif i%2==0 and not hold[i]: 
			callFitName+='p['+str(counter)+']'
			counter+=1
		else:
			callFitName+='*np.exp(-timepoints/p['+str(counter)+'])'
			counter+=1
		if i<len(hold)-1 and i%2==1:
			callFitName+='+'
	if bounds==0:
		callFitName+=')*np.heaviside(timepoints,0), timepoints, pixelSeries, p0=p1)'
	else:
		callFitName+=')*np.heaviside(timepoints,0), timepoints, pixelSeries, p0=p1,bounds=bounds)'
	return callFitName

def fitPixel(timepoints,pixelSeries,p0,p1,bounds,hold,callFitName):
	
	popt, pcov = eval(callFitName)
	for i in range(0, len(hold)):
		if hold[i]:
			popt=np.insert(popt,i,p0[i])
	yfit = ff.multiexp(timepoints, popt)
	return popt, yfit

def extractOscillations(timepoints,data,guessCoefs,bounds,hold):
	print 'we will be fitting to an exponential decay of order '+str(len(guessCoefs)/2)
	#define places to store the results
	start_time = time.time()
	callFitName=genCallFitName(guessCoefs,hold,bounds)
	holdloc=[i for i, x in enumerate(hold) if x]
	p1=np.delete(guessCoefs,holdloc)
	poptList=np.array([]);fitMat=np.array([]);resMat=np.array([]);
	#sweep through each pixel of the data fitting to an exponential and extract the residuals and organize results in a new matrix
	for i in range(0,data.shape[0]):
		#popt,yfit = fitPixel(timepoints,data[i,:],guessCoefs,bounds)
		popt,yfit = fitPixel(timepoints,data[i,:],guessCoefs,p1,bounds,hold,callFitName)
		poptList = np.append(poptList,popt)
		fitMat = np.append(fitMat,yfit)
		yres = data[i,:]-yfit
		resMat =np.append(resMat,yres)

	#lets keep our arrays organized.
	poptList=np.reshape(poptList,(data.shape[0],len(guessCoefs)))
	fitMat = np.reshape(fitMat,(data.shape[0],data.shape[1]))
	resMat = np.reshape(resMat,(data.shape[0],data.shape[1]))

	print "--- %s seconds --- to fit" % (time.time() - start_time),
	print str(data.shape[0])+" time traces"
	return poptList, fitMat, resMat

def removeBaseline(resMat_FFT,w_FSRS,polyorder):
	print 'we will be fitting the other FFT signal in the other dimension to polynomial of order '+str(polyorder)
	start_time = time.time()
	polyCoefsList=np.array([]);BaselineMat=np.array([]);ESFSRSMat=np.array([]);
	p0=np.array(np.ones(polyorder)*np.amax(resMat_FFT)/2.0)
	for i in range(0, resMat_FFT.shape[1]):
		popt, pcov = curve_fit(ff.polyn, w_FSRS, resMat_FFT[:,i], p0=p0)
		yfit = ff.polynomial(w_FSRS,popt)
		yres = resMat_FFT[:,i]-yfit
		polyCoefsList = np.append(polyCoefsList,popt)
		BaselineMat = np.append(BaselineMat,yfit)
		ESFSRSMat =np.append(ESFSRSMat,yres)
	polyCoefsList=np.reshape(polyCoefsList,(resMat_FFT.shape[1],len(p0)))
	BaselineMat = np.transpose(np.reshape(BaselineMat,(resMat_FFT.shape[1],resMat_FFT.shape[0])))
	ESFSRSMat = np.transpose(np.reshape(ESFSRSMat,(resMat_FFT.shape[1],resMat_FFT.shape[0])))
	print "--- %s seconds --- to fit" % (time.time() - start_time),
	print str(resMat_FFT.shape[1])+" frequency traces"
	return polyCoefsList, BaselineMat, ESFSRSMat