#python "C:\Users\skell\Documents\Indigo\TIn\Gaussian\DuschinskyAnalysis\duschinsky9.py"
#https://strawberryfields.ai/photonics/apps/run_tutorial_excitations.html
import numpy as np
import strawberryfields as sf
from strawberryfields.apps import qchem
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from itertools import islice
import os, sys
import pandas as pd
import math
#np.set_printoptions(threshold=sys.maxsize)

def findcom(r):
	com=r.mean(axis=0)
	return com

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def findlongaxis(r):
	a=r.mean(axis=0)
	uu, dd, vv = np.linalg.svd(r - a)
	# Now vv[0] contains the first principal component, i.e. the direction
	# vector of the 'best fit' line in the least squares sense.

	# Now generate some points along this best fit line, for plotting
	longaxis=np.vstack((a,a+vv[0,:]*10))
	return longaxis

def findlongvector(r):
	a=r.mean(axis=0)
	uu, dd, vv = np.linalg.svd(r - a)
	# Now vv[0] contains the first principal component, i.e. the direction
	# vector of the 'best fit' line in the least squares sense.
	print(vv)
	return np.vstack((vv,-vv))

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2)
    b, c, d = -axis*math.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def reflection_matrix(normunitvector):
    a=normunitvector[0]
    b=normunitvector[1]
    c=normunitvector[2]
    return np.array([[1-2*a**2, -2*a*b, -2*a*c],
                     [-2*a*b, 1-2*b**2, -2*b*c],
                     [-2*a*c, -2*b*c, 1-2*c**2]])


def getPerpendicularVector(v1,v2):
	pv=np.vstack((-np.cross(v1,v2),np.cross(v1,v2)))
	pv/=(pv[0,0]**2+pv[0,1]**2+pv[0,2]**2)**.5
	return pv


def buildMassDict():
	massDict={1:1.00782,
	6:12,
	8:15.99,
	16:32.065
	}
	return massDict
def getCoordLineIndex(fn,verbose=0):
	found=False
	start_Data_Index=-1
	end_Data_Index=-1
	with open(fn, 'r') as f:
		for index, line in enumerate(f):       
			# search string
			if not found:
				if ' Number     Number       Type             X           Y           Z' in line:
					if verbose:
						print('start string found in a file')
						print("index: ", index)
					start_Data_Index	=index+2
					found=True
					# don't look for next lines
			elif found:
				if '                    Distance matrix (angstroms):' in line:
					if verbose:
						print('end string found in a file')
						print("index: ", index)
					end_Data_Index	=index-1
					break
	print(start_Data_Index	,end_Data_Index)
	return start_Data_Index	,end_Data_Index

def getCoordData(fn,si,ei):
	d=[]
	with open(fn) as f:
	    for line in islice(f, si, ei):
	        d.append([float(x.replace('\n', '')) for x in line.split(" ") if not x==''])
	d_prelim=np.array(d)
	an_prelim=d_prelim[:,1]
	massDict=buildMassDict()
	m_prelim=[massDict[i] for i in an_prelim]
	massi=np.reshape(np.transpose(np.multiply(np.ones((3,len(m_prelim))),m_prelim)),(1,np.shape(m_prelim)[0]*3))[0]
	coordsi_prelim=d_prelim[:,3:6]
	coordsi=np.reshape(coordsi_prelim, (1,np.shape(coordsi_prelim)[0]*3))[0]
	return massi,coordsi,m_prelim,coordsi_prelim,an_prelim

def getFreqLineIndex(fn,verbose=0):
	freqIndexes=[]
	with open(fn, 'r') as f:
		for index, line in enumerate(f):       
			# search string
			if ' Frequencies --    ' in line:
				if verbose:
					print('start string found in a file')
					print("index: ", index)
					print(line)
				freqIndexes.append(index)
			elif ' Frequencies --   ' in line:
				if verbose:
					print('more lines?!')
					print("index: ", index)
					print(line)
				freqIndexes.append(index)
	return freqIndexes

def findbonds(coords,an_prelim,min_bl=2.1):
	bc=[]
	bl=[]
	for i in range(np.shape(coords)[0]):
		for j in range(i,np.shape(coords)[0]):
			if i != j:
				if an_prelim[i]!=1 and an_prelim[j]!=1:
					bondlength=((coords[i,0]-coords[j,0])**2+(coords[i,1]-coords[j,1])**2+(coords[i,2]-coords[j,2])**2)**.5
					if bondlength<min_bl:
						bl.append(bondlength)
						bc.append([i,j])
	print(bc)
	print(range(np.shape(coords)[0]))
	return bc,bl

def getFreqData(fn,fi):
	d=[]
	with open(fn) as f:
		for index, line in enumerate(f):
			if index in fi:
				d.append([float(x) for x in line.replace('Frequencies', '').replace('--', '').replace('\n','').split(" ") if not x==''])
				print(line)
		d2=np.array([d])
		d3=np.reshape(d2,np.shape(d2)[1]*3)
	return d3

def getNormalModes(fn,fi,natoms,datalineoffset,verbose=0):
	d=[]
	print([i for i in fi])
	fj=np.ndarray.flatten(np.array([np.arange(i+datalineoffset,i+datalineoffset+natoms) for i in fi]))
	with open(fn) as f:
		for index, line in enumerate(f):
			if index in fj:
				d.append([float(x.replace('\n', '')) for x in line.split(" ") if not x==''])
	d3=np.delete(d,[0,1],axis=1)
	d5= [[d3[j*natoms:(j+1)*natoms,3*i:3+3*i].T for i in range(3)] for j in range(26)]
	index=0
	d6=np.array([], dtype=np.float64).reshape(0,3)
	for i in range(np.shape(d5)[0]):
		for j in range(np.shape(d5)[1]):
			d6=np.vstack((d6,d5[i][j][:][:].T))
			if verbose:
				print(d5[i][j][:][:])
				print(d6)
	#print(d5[0][0][:][:])
	#print(d6[0:28,:])
	np.set_printoptions(threshold=sys.maxsize)
	#print(np.ravel(d6[0:28,:].T,order='F'))
	#print(np.reshape(np.ravel(d6[:,:].T,order='F'),(3*natoms,3*natoms-6)))
	d7=np.swapaxes(np.reshape(np.ravel(d6.T,order='F'),(3*natoms-6,3*natoms)),0,1)
	if verbose:
		print(d7)
		print(np.shape(d7))
	return d7

def extractGaussianVibrationData(fn,datalineoffset):
	start_Data_Index,end_Data_Index=getCoordLineIndex(fn)
	mass,coords,mass_prelim,coords_prelim,an_prelim=getCoordData(fn,start_Data_Index,end_Data_Index)
	natoms=np.shape(mass_prelim)[0]
	freqIndexes=getFreqLineIndex(fn,verbose=0)
	freq = getFreqData(fn,freqIndexes)
	print(freqIndexes)
	print("datalineoffset	:", datalineoffset)
	nmcoords = getNormalModes(fn,freqIndexes,natoms,datalineoffset,verbose=0)
	start_Data_Index=[]
	end_Data_Index=[]
	freqIndexes=[]
	return mass,coords,freq,nmcoords,mass_prelim,coords_prelim,an_prelim

natoms=28
fn_ex='C:\\Users\\skell\\Documents\\Indigo\\TIn\\Gaussian\\DuschinskyAnalysis\\Ex_Trans.log'
fn_gr='C:\\Users\\skell\\Documents\\Indigo\\TIn\\Gaussian\\DuschinskyAnalysis\\Gr_Trans.log'
destination="C:\\Users\\skell\\Documents\\Indigo\\TIn\\Gaussian\\DuschinskyAnalysis\\DeuchinskyMatrix_ttin_MK_Modes0.tif"
destination2="C:\\Users\\skell\\Documents\\Indigo\\TIn\\Gaussian\\DuschinskyAnalysis\\DeuchinskyMatrix_coords_ttin2.pdf"          
datalineoffset=5
mass_ex,coords_ex,freq_ex,nmcoords_ex,mass_prelim_ex,coords_prelim_ex,an_prelim_ex=extractGaussianVibrationData(fn_ex,datalineoffset)
datalineoffset=8
mass_gr,coords_gr,freq_gr,nmcoords_gr,mass_prelim_gr,coords_prelim_gr,an_prelim_gr=extractGaussianVibrationData(fn_gr,datalineoffset)
freq_gr_str=[format(q, '.0f') for q in freq_gr]
freq_ex_str=[format(r, '.0f') for r in freq_ex]
verbose=0
#reorient excited state molecule such that central c-c bond lies along y
nonhydrogen_index=[i for i in range(len(mass_prelim_gr)) if mass_prelim_gr[i]!= 1.00782]
carbon_index=[i for i in range(len(an_prelim_gr)) if an_prelim_gr[i]== 6]
oxygen_index=[i for i in range(len(an_prelim_gr)) if an_prelim_gr[i]== 8]
sulfur_index=[i for i in range(len(an_prelim_gr)) if an_prelim_gr[i]== 16]
centercarbons_index=[7,12]
carbon_bromine_bond_index=[[6,27],[13,26]]
coords_prelim_ex2=np.subtract(coords_prelim_ex,np.mean(coords_prelim_ex[[7,12],:],axis=0))
coords_prelim_gr2=np.subtract(coords_prelim_gr,np.mean(coords_prelim_gr[[7,12],:],axis=0))
acc=angle(coords_prelim_ex2[7,:], coords_prelim_gr2[7,:])
pv_acc=getPerpendicularVector(coords_prelim_ex2[7,:], coords_prelim_gr2[7,:])
abc=angle(coords_prelim_ex2[27,:], coords_prelim_gr2[27,:])
coords_prelim_ex_rot=np.multiply(np.dot(rotation_matrix(pv_acc[0,:],-acc), coords_prelim_ex2.T).T,1)
coords_prelim_ex_rot=np.multiply(np.dot(rotation_matrix(coords_prelim_ex_rot[7,:],abc-19.24*np.pi/180), coords_prelim_ex_rot.T).T,1)
coords_ex=np.reshape(coords_prelim_ex2, (1,np.shape(coords_prelim_ex2)[0]*3))[0]
coords_gr=np.reshape(coords_prelim_gr2, (1,np.shape(coords_prelim_gr2)[0]*3))[0]
coords_gr_rot=np.reshape(coords_prelim_gr2, (1,np.shape(coords_prelim_gr2)[0]*3))[0]
coords_ex_rot=np.reshape(coords_prelim_ex_rot, (1,np.shape(coords_prelim_ex_rot)[0]*3))[0]
bonds_ex,bl_ex=findbonds(coords_prelim_ex_rot,an_prelim_ex)
bonds_gr,bl_gr=findbonds(coords_prelim_gr2,an_prelim_gr)
delta_bl=np.subtract(bl_ex,bl_gr)
bondlabeloffset=[[0.25,-0.05],[-.1,-.05],[0.0,0],[0.1,0],[-0.25,0],[0.3,-0.1],[0.2,0],[0.65,-.2],[0.15,0],[0.2,0],[0.1,.2],[0.15,0],[0.05,-.05],[0.1,0.3],[0.2,-0.1],[0.1,0],[0.2,0],[0.2,0],[.3,.1],[0.1,-.1],[-.1,.1],[-.05,.05],[0.15,-.25]]
#                     0           1        2       3       4       5          6         7         8        9       10        11         12        13         14       15     16      17       18       19      20       21        22 
print("sample:",coords_prelim_ex_rot[bonds_ex[6]])
midbondcoord_ex=np.vstack([np.mean(coords_prelim_ex_rot[bonds_ex[i]],axis=0) for i in range(np.shape(bonds_ex)[0])])
midbondcoord_gr=np.vstack([np.mean(coords_prelim_gr2[bonds_gr[i]],axis=0) for i in range(np.shape(bonds_gr)[0])])
midbondcoord_ex2=np.reshape(midbondcoord_ex, (1,np.shape(midbondcoord_ex)[0]*3))[0]
print("bl",bl_ex)
print("middle bond coords",midbondcoord_ex)
xdatai=coords_gr_rot[0::3]
ydatai=coords_gr_rot[1::3]
zdatai=coords_gr_rot[2::3]
xdataf=coords_ex_rot[0::3]
ydataf=coords_ex_rot[1::3]
zdataf=coords_ex_rot[2::3]
if 1:
	Ud, delta = qchem.duschinsky(nmcoords_gr, nmcoords_ex, coords_gr_rot, coords_ex_rot, freq_ex, mass_gr)
	Ud2=np.where(Ud > .4, Ud, 0)**1.5
	ax = plt.axes()
	plt.imshow(abs(Ud2), cmap="Greys")
	plt.colorbar(pad=0.15)
	plt.xlabel("Mode index")
	plt.ylabel("Mode index")
	filled_marker_style8 = dict(marker='None', linestyle='solid', markersize=5,color='red',markerfacecolor='tab:blue',markerfacecoloralt='lightsteelblue',markeredgecolor='red')
	ax.annotate("Frequency S1 (cm"+r'$^{-1}'+')', xy=(1.17, .5), xycoords='axes fraction', xytext=(1.17,.5),arrowprops=dict(arrowstyle="-", color='k'),rotation=90, ha='left', va='center')
	ax.annotate("Frequency S1 (cm"+r'$^{-1}'+')', xy=(.5, 1.17), xycoords='axes fraction', xytext=(.5, 1.17),arrowprops=dict(arrowstyle="-", color='k'),rotation=0, ha='center', va='bottom')
	for i in range(np.shape(Ud2)[0])[0::6]:
		ax.annotate(freq_gr_str[i], xy=((i+1)/np.shape(Ud2)[0], 1.0), xycoords='axes fraction', xytext=((i+1)/np.shape(Ud2)[0], 1.04),arrowprops=dict(arrowstyle="-", color='k'),rotation=90, ha='center', va='bottom')
		#ax.annotate(freq_gr_str[i], ((i-1)/np.shape(Ud2)[0], 1.07), xycoords='axes fraction', rotation=90)
	for j in range(np.shape(Ud2)[1])[0::6]:
		ax.annotate(freq_ex_str[j], xy=(1.0,1-(j+1)/np.shape(Ud2)[1]), xycoords='axes fraction', xytext=(1.04,1-(j+1)/np.shape(Ud2)[1]),arrowprops=dict(arrowstyle="-", color='k'), ha='left', va='center')
		#ax.annotate(freq_ex_str[i], (1.07,(j-1)/np.shape(Ud2)[1]), xycoords='axes fraction', rotation=0)
	ax.set_ylim([np.shape(Ud2)[1],0])
	plt.tight_layout()
	plt.savefig(destination)
	plt.show()
if 0:
	Ud, delta = qchem.duschinsky(nmcoords_gr, nmcoords_ex, coords_gr_rot, coords_ex_rot, freq_ex, mass_gr)
	Ud2=np.where(Ud > .15, Ud, 0)**100
	fig, ax = plt.subplots(constrained_layout=True)
	#ax.get_xaxis().set_visible(False)
	#ax1 = ax.twiny()
	ax2 = ax.twinx().twiny()
	#ax2.plot(x, y2, 'b-')
	#ax2.xaxis.set_label_position('bottom') 
	#ax1.set_xlabel('TOP', labelpad=20)
	#ax2.set_xlabel('BOTTOM', labelpad=20)
	# Set equal limits on both yaxis so that the ticks line up
	#ax.set_ylim(ax1.get_ylim())
	# Set the tick locations and labels
	ax2.set_xlim([1, np.shape(Ud2)[0]])
	ax2.set_ylim([np.shape(Ud2)[1],1])
	ax.set_xticks(range(np.shape(Ud2)[0])[0::6],labels=freq_gr_str[0::6])
	ax.set_yticks(range(np.shape(Ud2)[1])[0::6],labels=freq_ex_str[0::6])	
	ax.imshow(abs(Ud2), cmap="Greens")
	#ax2.imshow(abs(Ud2), cmap="Greens")
	#ax3.imshow(abs(Ud2), cmap="Greens")
	#pcm = ax.pcolormesh(range(np.shape(Ud)[0]),range(np.shape(Ud)[1]),Ud)
	#fig.colorbar(pcm)
	fig.colorbar(mappable=None,cmap='Greens',ax=ax, orientation='vertical', label='a colorbar label')
	#plt.xlabel("Mode index")
	#plt.ylabel("Mode index")
	#plt.tight_layout()
	#ax1.set_xticks(range(np.shape(Ud)[0])[0::6])
	#ax1.set_xticklabels(freq_gr_str[0::6])
	#ax1.set_yticks(range(np.shape(Ud)[1])[0::6])
	#ax1.set_yticklabels(freq_ex_str[0::6])
	plt.savefig(destination)
	plt.show()
	# Data for three-dimensional scattered points
if 0:
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.scatter3D(xdatai[nonhydrogen_index], ydatai[nonhydrogen_index], zdatai[nonhydrogen_index], color='k');
	ax.scatter3D(xdataf[nonhydrogen_index], ydataf[nonhydrogen_index], zdataf[nonhydrogen_index], color='red');
	ax.set_xlim3d(-6,6)
	ax.set_ylim3d(-6,6)
	ax.set_zlim3d(-6,6)
	ax.set_xlabel('$X$', fontsize=10)
	ax.set_ylabel('$Y$', fontsize=10)
	plt.show()
if verbose:
	print("nmcoords_gr:",np.shape(nmcoords_gr))
	print("nmcoords_ex:",np.shape(nmcoords_ex))
	print("coords_gr:",np.shape(coords_gr))
	print("coords_ex:",np.shape(coords_ex))
	print("freq:",np.shape(freq_ex))
	print("mass:",np.shape(mass_gr))
	print("mass:",mass_prelim_gr)
	print("coords_prelim_ex",coords_prelim_ex2[[7,12]])
	print("coords_prelim_gr",coords_prelim_gr2[[7,12]])
	print("acc: ", acc)
	print("abc: ",abc)
	print("coords_prelim_ex_rot",coords_prelim_ex_rot)
	print("coords_prelim_gr_rot",coords_prelim_gr2)
if 0:
	#fig = plt.figure(figsize=(6.66,2.22))
	fig = plt.figure()
	#fig.suptitle('thioindigio excited state coordnates', fontsize=14)
	#fig.subplots_adjust(left=0.4)
	filled_marker_style0 = dict(marker='o', linestyle='None', markersize=12,color='black',markerfacecolor='white',markerfacecoloralt='white',markeredgecolor='black')
	filled_marker_style1 = dict(marker='$C$', linestyle='None', markersize=5,color='darkgrey',markerfacecolor='white',markerfacecoloralt='white',markeredgecolor='black')
	filled_marker_style2 = dict(marker='o', linestyle='None', markersize=12,color='red',markerfacecolor='white',markerfacecoloralt='white',markeredgecolor='red')
	filled_marker_style3 = dict(marker='$C$', linestyle='None', markersize=5,color='darkgrey',markerfacecolor='tab:blue',markerfacecoloralt='lightsteelblue',markeredgecolor='red')
	filled_marker_style4 = dict(marker='$O$', linestyle='None', markersize=5,color='darkgrey',markerfacecolor='tab:blue',markerfacecoloralt='lightsteelblue',markeredgecolor='black')
	filled_marker_style5 = dict(marker='$O$', linestyle='None', markersize=5,color='darkgrey',markerfacecolor='tab:blue',markerfacecoloralt='lightsteelblue',markeredgecolor='red')
	filled_marker_style6 = dict(marker='$S$', linestyle='None', markersize=5,color='darkgrey',markerfacecolor='tab:blue',markerfacecoloralt='lightsteelblue',markeredgecolor='black')
	filled_marker_style7 = dict(marker='$S$', linestyle='None', markersize=5,color='darkgrey',markerfacecolor='tab:blue',markerfacecoloralt='lightsteelblue',markeredgecolor='red')
	filled_marker_style8 = dict(marker='None', linestyle='dashed', markersize=5,color='red',markerfacecolor='tab:blue',markerfacecoloralt='lightsteelblue',markeredgecolor='red')
	filled_marker_style9= dict(marker='None', linestyle='solid', markersize=5,color='black',markerfacecolor='tab:blue',markerfacecoloralt='lightsteelblue',markeredgecolor='black')
	filled_marker_style10= dict(marker='None', linestyle='solid', markersize=5,color='black',markerfacecolor='tab:blue',markerfacecoloralt='lightsteelblue',markeredgecolor='black')
	ax = plt.axes()
	plt.rcParams["font.family"] = "serif"
	plt.rcParams["font.serif"] = ["Times New Roman"]
	for i in range(np.shape(bonds_ex)[0]):
		print(i)
		#print(midbondcoord_gr[i,1],bondlabeloffset),np.add(midbondcoord_gr[i,0],bondlabeloffset[i,1])
		ax.plot(ydataf[bonds_ex[i]],xdataf[bonds_ex[i]], fillstyle='full', **filled_marker_style8);
		ax.plot(ydatai[bonds_gr[i]],xdatai[bonds_gr[i]], fillstyle='full', **filled_marker_style9);
		#currentlabel=str(i)+" "+format(bl_gr[i], '.2f')+"  ("+format(delta_bl[i], '.2f')+")"
		currentlabel=format(bl_gr[i], '.2f')+"  ("+format(delta_bl[i], '.2f')+")"
		ax.text(np.add(midbondcoord_gr[i,1],bondlabeloffset[i][0]),np.add(midbondcoord_gr[i,0],bondlabeloffset[i][1]), currentlabel, fontsize=10,color='black', weight='bold',verticalalignment='center', ha='center')
		#ax.text(midbondcoord_gr[i,1]+bondlabeloffset_x_gr[i],midbondcoord_gr[i,0]+bondlabeloffset_y_gr[i], currentlabel, fontsize=8,color='black', weight='bold',verticalalignment='center', ha='center')
		#ax.text(midbondcoord_ex[i,1]+bondlabeloffset_x_ex[i],midbondcoord_ex[i,0]+bondlabeloffset_y_ex[i], str(i)+"("+format(bl_ex[i], '.2f')+")", fontsize=8,color='red', weight='bold',verticalalignment='center', ha='center')
	ax.plot(ydatai[nonhydrogen_index],xdatai[nonhydrogen_index], fillstyle='full', **filled_marker_style0);
	ax.plot(ydataf[nonhydrogen_index],xdataf[nonhydrogen_index], fillstyle='full', **filled_marker_style2);
	ax.plot(ydatai[carbon_index],xdatai[carbon_index], fillstyle='full', **filled_marker_style1);
	ax.plot(ydataf[carbon_index],xdataf[carbon_index], fillstyle='full', **filled_marker_style3);
	ax.plot(ydatai[oxygen_index],xdatai[oxygen_index], fillstyle='full', **filled_marker_style4);
	ax.plot(ydataf[oxygen_index],xdataf[oxygen_index], fillstyle='full', **filled_marker_style5);
	ax.plot(ydatai[sulfur_index],xdatai[sulfur_index], fillstyle='full', **filled_marker_style6);
	ax.plot(ydataf[sulfur_index],xdataf[sulfur_index], fillstyle='full', **filled_marker_style7);
	ax.plot([4,5],[-2,-2], fillstyle='full', linewidth=3, **filled_marker_style10);
	ax.text(4.5,-2.15, "1 nm", fontsize=11,color='black', weight='bold',verticalalignment='top', ha='center')
	fig.tight_layout()
	plt.axis('scaled')
	plt.axis('off')
	plt.savefig(destination2, format="pdf", bbox_inches="tight")
	#plt.savefig(destination2)
	plt.show()
