import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import linalg
from numpy import matrix
import matplotlib.pyplot as pl
import init as cc
mpl.rcParams['legend.fontsize'] = 10


#use mercury to export .xyz file to carbonfile with packing turned on
#angles and lengths of unit cell atomic coordiantes file name
alpha=(101.245)/180*(np.pi)
beta=(99.432)/180*(np.pi)
gamma=(94.272)/180*(np.pi)

a=6.0394
b =7.8075
c =12.528


#csc(x)=1/sin(x)   ,   cot=1/tan(x)
p=np.zeros((8,3))
#define origin and points of unitcell
p[0,0]=p[0,1]=p[0,2]=0
p[1,0]=p[0,0]+c*np.cos(beta)
p[1,1]=p[0,1]+c*(np.cos(gamma)*np.sin(alpha)**-1-np.cos(beta)*np.tan(alpha)**-1)
p[1,2]=p[0,2]+c*(np.sin(beta)**2+(np.cos(gamma)*np.cos(alpha)**-1-np.cos(beta)*np.tan(alpha)**-1)**2)**.5
p[2,0]=p[0,0]+b*np.cos(alpha)
p[2,1]=p[0,1]+b*np.sin(alpha)
p[2,2]=p[0,2]
p[3,0]=p[0,0]+a
p[3,1]=p[0,1]
p[3,2]=p[0,2]
p[4,0]=p[1,0]+p[2,0]
p[4,1]=p[1,1]+p[2,1]
p[4,2]=p[1,2]+p[2,2]
p[5,0]=p[1,0]+p[3,0]
p[5,1]=p[1,1]+p[3,1]
p[5,2]=p[1,2]+p[3,2]
p[6,0]=p[2,0]+p[3,0]
p[6,1]=p[2,1]+p[3,1]
p[6,2]=p[2,2]+p[3,2]
p[7,0]=p[1,0]+p[2,0]+p[3,0]
p[7,1]=p[1,1]+p[2,1]+p[3,1]
p[7,2]=p[1,2]+p[2,2]+p[3,2]
mid=p[7,0:]/2
#Load Carbon and hydrogen coordinates
carbonfile='tccoordinates.txt'
hydrogenfile='hydrogenccoordinates.txt'
coords=np.genfromtxt(carbonfile,delimiter="	")
hcoords=np.genfromtxt('hydrogencoordinates.txt',delimiter="	")
call=cc.findcom(coords[0:,:])
offset=mid-call


#recenter unit cell around middle of unitcell
p2=np.subtract(p,mid)
#center actomic coordinates around center of mass/point of inversion:
coords2=np.subtract(coords,call)


hcoords2=np.subtract(hcoords,call)
l1=np.subtract(cc.findlongaxis(coords[0:18,:]),call)
l2=np.subtract(cc.findlongaxis(coords[18:36,:]),call)
l3=np.subtract(cc.findlongaxis(coords[36:54,:]),call)
l4=np.subtract(cc.findlongaxis(coords[54:72,:]),call)
v0=cc.findlongvector(coords2[0:18,:])
v1=cc.findlongvector(coords2[18:36,:])
c1=np.subtract(cc.findcom(coords[0:18,:]),call)
c2=np.subtract(cc.findcom(coords[18:36,:]),call)
c3=np.subtract(cc.findcom(coords[36:54,:]),call)
c4=np.subtract(cc.findcom(coords[54:72,:]),call)
com=np.vstack((c1,c2,c3,c4))
mid2=[0,0,0]
coords3=coords2.reshape(4,18,3)


fig = pl.figure()
ax = fig.gca(projection='3d')
cinv=np.multiply(coords2,-1)
theta=np.pi
crot=np.multiply(np.dot(cc.rotation_matrix(v1[0,:],theta), coords2.T).T,1)
cref=np.multiply(np.dot(cc.reflection_matrix(v1[0,:]), coords2.T).T,1)
cs2=np.dot(cc.reflection_matrix(v1[0,:]), np.dot(cc.rotation_matrix(v1[0,:],theta), coords2.T)).T
#ax.scatter(coords2[0:,0],coords2[0:,1],coords2[0:,2])
#ax.scatter(cinv[0:,0],cinv[0:,1],cinv[0:,2],color='m')
#ax.scatter(crot[0:,0],crot[0:,1],crot[0:,2],color='y')
ax.scatter(cref[0:,0],cref[0:,1],cref[0:,2],color='y')
#ax.scatter(mid[0],mid[1],mid[2],label='M',color='g')
ax.scatter(cs2[0:,0],cs2[0:,1],cs2[0:,2],color='c')



if 0:
	ax.plot(l1[:,0],l1[:,1],l1[:,2],color='r')
	ax.plot(l2[:,0],l2[:,1],l2[:,2],color='r')
	ax.plot(l3[:,0],l3[:,1],l3[:,2],color='r')
	ax.plot(l4[:,0],l4[:,1],l4[:,2],color='r')


if 1:
	ax.plot(10*v0[:,0],10*v0[:,1],10*v0[:,2],color='m')



if 0:
	ax.scatter(com[:,0],com[:,1],com[:,2],label='H',color='g')
	ax.scatter(mid2[0],mid2[1],mid2[2],label='J',color='m')


if 0:
	ax.scatter(hcoords2[0:,0],hcoords2[0:,1],hcoords2[0:,2],label='H',color='r')
if 1:
	#parallelogram 0
	ax.plot([p2[0,0],p2[1,0],p2[5,0],p2[3,0],p2[0,0]],[p2[0,1],p2[1,1],p2[5,1],p2[3,1],p2[0,1]],[p2[0,2],p2[1,2],p2[5,2],p2[3,2],p2[0,2]],color='b')
	#parallelogram 1
	ax.plot([p2[0,0],p2[1,0],p2[4,0],p2[2,0],p2[0,0]],[p2[0,1],p2[1,1],p2[4,1],p2[2,1],p2[0,1]],[p2[0,2],p2[1,2],p2[4,2],p2[2,2],p2[0,2]],color='b')
	#parallelogram 2
	ax.plot([p2[0,0],p2[2,0],p2[6,0],p2[3,0],p2[0,0]],[p2[0,1],p2[2,1],p2[6,1],p2[3,1],p2[0,1]],[p2[0,2],p2[2,2],p2[6,2],p2[3,2],p2[0,2]],color='b')
	#parallelogram 3
	ax.plot([p2[1,0],p2[4,0],p2[7,0],p2[5,0],p2[1,0]],[p2[1,1],p2[4,1],p2[7,1],p2[5,1],p2[1,1]],[p2[1,2],p2[4,2],p2[7,2],p2[5,2],p2[1,2]],color='b')
	#parallelogram 4
	ax.plot([p2[3,0],p2[5,0],p2[7,0],p2[6,0],p2[3,0]],[p2[3,1],p2[5,1],p2[7,1],p2[6,1],p2[3,1]],[p2[3,2],p2[5,2],p2[7,2],p2[6,2],p2[3,2]],color='b')
	#parallelogram 5
	ax.plot([p2[2,0],p2[4,0],p2[7,0],p2[6,0],p2[2,0]],[p2[2,1],p2[4,1],p2[7,1],p2[6,1],p2[2,1]],[p2[2,2],p2[4,2],p2[7,2],p2[6,2],p2[2,2]],color='b')

#add connectivity
if 1:
	for i in np.arange(4):
		ax.plot(coords3[i,0:,0],coords3[i,0:,1],coords3[i,0:,2],color='k')
		
		for j in np.arange(4):
			ax.plot(coords3[i,np.array([-1-2*j,2*j]),0],coords3[i,np.array([-1-2*j,2*j]),1],coords3[i,np.array([-1-2*j,2*j]),2],color='k')


#ax.legend()
ax.set_xlim3d(-7.5,7.5)
ax.set_ylim3d(-7.5,7.5)
ax.set_zlim3d(-7.5,7.5)

plotting = True
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

pl.show()

