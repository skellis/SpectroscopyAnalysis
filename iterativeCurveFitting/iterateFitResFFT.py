import init as cc
import graphing as gr
import fitfuncs as ff
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

if 1:
	#Section to load existing data
	dataname ='RGTCNQTMB.txt'
	timepointsName='timepointsTCNQTMB.txt'
	shiftxName='shiftxTCNQTMB.txt'
	data=np.genfromtxt(dataname,delimiter='	')
	timepoints=np.genfromtxt(timepointsName,delimiter='	')
	if 1:
		st=120 # fs
		data,timepoints = cc.trimData(data,timepoints,st)
	shiftx=np.genfromtxt(shiftxName,delimiter='	')
	t,w=np.meshgrid(timepoints,shiftx)

	print 'loading data from text files:',
	print dataname+', ',
	print timepointsName+', ',
	print shiftxName+', '
else:
	#Section to generate Dummy Data
	#If you don't want this component included in your data you can pass it zero amplitude or a null array. For example expParams=[]
	#note sometimes you get a runtimewarning when you inclue negative time points in fitting.
	timepoints = np.arange(0.,4000.,20.)
	shiftx=np.arange(250.,1250.,10.)
	#FSRSCoords=np.array([amp (mOD),phase (radians),freq (cm-1),fwhm (cm-1),tauPeak (fs)...])
	FSRSCoords=[1.,0/2,1000.,30.,2000.,2.,0/2,500.,30.,2000.]
	#ISRSCoords=np.array([amp,freq (cm-1),taudephasing (fs)])
	ISRSCoords=[1.,333.,1000.]
	#ExpParams=np.array([amp1,tau1...])
	expParams=[-3,200,6.,2000.] 
	noise=[.1]
	data,t,w=cc.generateDummyData(shiftx,timepoints,FSRSCoords,ISRSCoords,expParams,noise)
	if 0:
		np.savetxt('simtimepoints.txt',timepoints,delimiter='	',header='fs')
		np.savetxt('simshiftx.txt',shiftx,delimiter='	',header='cm-1')
		np.savetxt('simData.txt',data,delimiter='	',header='mOD')
	print 'Generating data from input parameters.'
#Begin the fitting
guessCoefs=[-0.02, 163.76,-0.03, 1259.24,.05,10300]
#guessCoefs=[-1,200,8.,3000.,2,4000]
#Hold must be of same length as guessCoefss.
#a zero indicates that the coefficient is to be fit while a 1 indicates that this parameter is held fixed to guessCoefs[i]
hold=[0,1,0,1,0,1]
#set the bounds for curve fitting.
#a good first guess for the amplitudes is probably 
#formatBounds(guessCoefs,minAmp,maxAmp,mintau,maxTau)
bounds=cc.formatBounds(guessCoefs,hold,-10*abs(np.amin(data)),10*abs(np.amax(data)),0,100000)
#minamp=-10.; maxamp=10.; mintimeconstant=0.; maxtimeconstant=100000.;
poptList,fitMat,resMat= cc.extractOscillations(timepoints,data,guessCoefs,bounds,hold)
pad=300
power=1
resMat_FFT=np.power(np.real(np.fft.fft(resMat,n=np.shape(resMat)[1]+pad,axis=1)),power)
#np.fft.fft yields a matrix of redundant data where the frequecy axis goes up to the nyquist limit and then back to zero.
#extract the portion of the FFT matrix we are interested in.
resMat_FFT = resMat_FFT[:,0:np.shape(resMat_FFT)[1]/2][:]
freqpoints=np.linspace(0.,cc.wn2fs(timepoints[1]-timepoints[0])/2,num=np.shape(resMat_FFT)[1])
wISRS,wFSRS=np.meshgrid(freqpoints,shiftx)
polyorder=6
polyCoefsList,BaselineMat,ESFSRSMat=cc.removeBaseline(resMat_FFT,shiftx,polyorder)


#Turn this on if would like to save the results
if 1:
	np.savetxt('2DESFSRSTCNQTMB_withBL.txt',resMat_FFT,delimiter='	',header='Spectral Density')
	np.savetxt('ResidualTCNQTMB.txt',resMat,delimiter='	',header='Spectral Density')
	np.savetxt('FitTCNQTMB.txt',fitMat,delimiter='	',header='Spectral Density')
	np.savetxt('ImpulsiveFreq.txt',freqpoints,delimiter='	',header='cm-1')
	np.savetxt('2DESFSRSTCNQTMB_noBL.txt',ESFSRSMat,delimiter='	',header='Spectral Density')

#Below are various plotting functions which will help you trouble shoot your data analysis.
if 1:
	#Functions for plotting a slice along the time domain from the FSRS data matrix
	desiredTime=500 #fs
	#plotTimeSlice(data,timepoints,shiftx,desiredTime,saveFig,autoClose)
	gr.plotTimeSlice(data,fitMat,timepoints,shiftx,desiredTime,0,0)
if 1:
	#A function for plotting a slice along the freq domain from the FSRS data matrix
	desiredFreq=882 #cm-1
	#plotFreqSlice(data,timepoints,shiftx,desiredFreq,saveFig,autoClose)
	gr.plotFreqSlice(data,fitMat,timepoints,shiftx,desiredFreq,0,0)
if 1:
	#plotFSRSContour(data,t,w,saveFig,autoClose)
	gr.plotFSRSContour(data,t,w,0,0)
if 0:
	#plotFSRSContour(fitMat,t,w,saveFig,autoClose)
	gr.plotFSRSContour(fitMat,t,w,0,0)
if 1:
	#plotFSRSContour(resMat,t,w,saveFig,autoClose)
	gr.plotFSRSContour(resMat,t,w,0,0)
if 1:
	#plotFSRSContour(resMat,t,w,saveFig,autoClose,minamp,maxamp)
	gr.plot2DContour(resMat_FFT.real,wISRS,wFSRS,0,0,-.02,.02)
if 1:
	gr.plotFitResults(poptList,timepoints,shiftx,0,0)
if 1:
	for i in range(0,int(np.floor(resMat_FFT.shape[1]/100))):
		plt.plot(shiftx, resMat_FFT[:,100*(i-1)],shiftx,BaselineMat[:,100*(i-1)])
	plt.show()
if 0:
	#plotFSRSContour(resMat,t,w,saveFig,autoClose,minamp,maxamp)
	gr.plot2DContour(BaselineMat,wISRS,wFSRS,0,0,-.02,.02)
if 1:
	#plotFSRSContour(resMat,t,w,saveFig,autoClose,minamp,maxamp)
	gr.plot2DContour(ESFSRSMat.real,wISRS,wFSRS,0,0,-.0025,.0025)