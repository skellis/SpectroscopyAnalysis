import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

def plotFitResults(poptList,timepoints,shiftx,saveFig,autoClose):

	plt.figure(1)
	for i in range(0,poptList.shape[1]/2):	
		ax=plt.subplot(poptList.shape[1]/2,2,2*i+1)
		ax.plot(shiftx, poptList[:,2*i])
		ax.set_xlabel('Raman Shift (cm-1)')
		ax.set_ylabel('Amplitude')
		ax.text(0.5, 0.9, 'Amplitude Component '+str(i+1), horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
		
		ax.text(0.5, 0.2, 'Median = '+str(np.around(np.median(poptList[:,2*i]),decimals=2)), horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
		
		
		ax=plt.subplot(poptList.shape[1]/2,2,2*i+2)
		ax.plot(shiftx, poptList[:,2*i+1])
		ax.set_xlabel('Raman Shift (cm-1)')
		ax.set_ylabel('Time Constant (fs)')
		ax.text(0.5, 0.9, 'Time Constant '+str(i+1), horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
		
		ax.text(0.5, 0.2, 'Median = '+str(np.around(np.median(poptList[:,2*i+1]),decimals=2)), horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
		ax.yaxis.tick_right()		
		
	if saveFig:
		plt.savefig(time.strftime("%Y%m%d")+'ampFreqFitCoefs.pdf')

	if autoClose:
		plt.show(block=False)
		plt.pause(3)
		plt.close()
	else:
		plt.show()

def plotFSRSContour(data,t,w,saveFig,autoClose):

	matplotlib.rcParams['xtick.direction'] = 'out'
	matplotlib.rcParams['ytick.direction'] = 'out'
	plt.figure()
	im = plt.imshow(data, interpolation='bilinear', origin='lower', cmap=cm.jet,extent=(t[0,0]	, t[-1,-1], w[0,0], w[-1,-1]))
	levels = np.arange(np.amin(data),np.amax(data),(np.amax(data)-np.amin(data))/4)

	CS = plt.contour(data, levels, origin='lower', linewidths=1, extent=(t[0,0]	, t[-1,-1], w[0,0], w[-1,-1]))
 	
	# Thicken the zero contour.
	zc = CS.collections
	plt.setp(zc, linewidth=1)

	plt.title('FSRS Time Evolution')
	plt.flag()

	# We can still add a colorbar for the image, too.
	CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)

	l, b, width, height = plt.gca().get_position().bounds

	if saveFig:
		plt.savefig('ContourFSRSEvolution0.pdf')
	
	
	if autoClose:
		plt.show(block=False)
		plt.pause(3)
		plt.close()
	else:
		plt.show()

def plot2DContour(data,t,w,saveFig,autoClose,vmin=0,vmax=0):
	if vmin==0:
		vmin=-abs(np.amin(data))
	if vmax==0:
		vmax=abs(np.amax(data))
	else:
		print 'If the image plot is not mostly white it is recommended that you change vmin,vmax to defaults'
	matplotlib.rcParams['xtick.direction'] = 'out'
	matplotlib.rcParams['ytick.direction'] = 'out'
	plt.figure()
	im = plt.imshow(data, interpolation='bilinear', origin='lower', cmap=cm.bwr, vmin=vmin, vmax=vmax,extent=(t[0,0]	, t[-1,-1], w[0,0], w[-1,-1]))
	levels = np.arange(np.amin(data),np.amax(data),(np.amax(data)-np.amin(data))/4)
	levels=[]
	#CS = plt.contour(data, levels, origin='lower', linewidths=1, extent=(t[0,0]	, t[-1,-1], w[0,0], w[-1,-1]))
 	
	# Thicken the zero contour.
	#zc = CS.collections
	#plt.setp(zc, linewidth=1)

	#plt.title('FSRS Time Evolution')
	#plt.flag()

	# We can still add a colorbar for the image, too.
	#CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)

	#l, b, width, height = plt.gca().get_position().bounds

	if saveFig:
		plt.savefig('ContourFSRSEvolution0.pdf')
	
	
	if autoClose:
		plt.show(block=False)
		plt.pause(3)
		plt.close()
	else:
		plt.show()


def plotFreqSlice(data,fitMat,timepoints,shiftx,desiredFreq,saveFig,autoClose):
	
	if desiredFreq%(shiftx[1]-shiftx[0])!=0:
		print 'The desired frequency doesnt fall on a measured frequency so the previous frequency point will be displayed.'

	desiredLoc=(int(desiredFreq)-int(shiftx[0]))/(int(shiftx[1])-int(shiftx[0]))
	slice1=data[desiredLoc,:]
	fitslice1=fitMat[desiredLoc,:]
	matplotlib.rcParams['axes.unicode_minus'] = False
	fig, ax = plt.subplots()
	ax.plot(timepoints, slice1,timepoints,fitslice1)
	ax.set_title('FSRS signal at '+str(desiredFreq)+' cm^-1')
	if saveFig:
		plt.savefig('slice'+str(desiredFreq)+'cm.pdf')

	if autoClose:
		plt.show(block=False)
		plt.pause(3)
		plt.close()
	else:
		plt.show()


def plotTimeSlice(data,fitMat,timepoints,shiftx,desiredTime,saveFig,autoClose):

	if desiredTime%(timepoints[1]-timepoints[0])!=0:
		print 'The desired time doesnt fall on a measured timepoint so the previous timepoint point will be displayed.'
	desiredLoc=(int(desiredTime)-int(timepoints[0]))/(int(timepoints[1])-int(timepoints[0]))
	slice0=data[:,desiredLoc]
	fitslice0=fitMat[:,desiredLoc]
	matplotlib.rcParams['axes.unicode_minus'] = False
	fig, ax = plt.subplots()
	ax.plot(shiftx, slice0,shiftx,fitslice0)
	ax.set_title('FSRS signal at '+str(desiredTime)+' fs')
	if saveFig:
		plt.savefig('slice'+str(desiredTime)+'fs.pdf')
	if autoClose:
		plt.show(block=False)
		plt.pause(3)
		plt.close()
	else:
		plt.show()