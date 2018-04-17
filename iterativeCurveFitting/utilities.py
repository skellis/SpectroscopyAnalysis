import numpy as np

def multiexp2(t,*args):
	val=np.zeros(len(t))
	for i in np.arange(0,len(args),2):
		val+=args[i]*np.exp(-t/args[i+1])*np.heaviside(t,0)
	return val

def multiexp(t,args):
	val=np.zeros(len(t))
	for i in np.arange(0,len(args),2):
		val+=args[i]*np.exp(-t/args[i+1])*np.heaviside(t,0)
	return val

def genCallFitName(p0,hold):
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
	callFitName+=')*np.heaviside(timepoints,0), timepoints, y, p0=p1)'
	return callFitName