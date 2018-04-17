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

def exp1(t,a1,tau1):
	val=0.
	val=(a1*np.exp(-t/tau1))*np.heaviside(t,0)
	return val

def exp2(t,a1,tau1,a2,tau2):
	val=0.
	val=(a1*np.exp(-t/tau1)+a2*np.exp(-t/tau2))*np.heaviside(t,0)
	return val

def exp3(t,a1,tau1,a2,tau2,a3,tau3):
	val=0.
	val=(a1*np.exp(-t/tau1)+a2*np.exp(-t/tau2)+a3*np.exp(-t/tau3))*np.heaviside(t,0)
	return val

def exp4(t,a1,tau1,a2,tau2,a3,tau3,a4,tau4):
	val=0.
	val=(a1*np.exp(-t/tau1)+a2*np.exp(-t/tau2)+a3*np.exp(-t/tau3)+a3*np.exp(-t/tau3))*np.heaviside(t,0)
	return val

def polyn(w,*p):
	val=np.zeros(len(w))
	for i in range(0, len(p)):
		val+=p[i]*w**i
	return val

def polynomial(w,coefs):
	val=np.zeros(len(w))
	for i in range(0, len(coefs)):
		val+=coefs[i]*w**i
	return val