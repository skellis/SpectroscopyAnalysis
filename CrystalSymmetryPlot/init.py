import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from scipy import linalg
from numpy import matrix
import matplotlib.pyplot as pl
import sys
import glob
import os


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

	# Now generate some points along this best fit line, for plotting.

	# I use -7, 7 since the spread of the data is roughly 14
	# and we want it to have mean 0 (like the points we did
	# the svd on). Also, it's a straight line, so we only need 2 points.
	linepts = vv[0] * np.mgrid[-7:7:2j][:, np.newaxis]
	
	# shift by the mean to get the line in the right place
	linepts += a

	# Verify that everything looks right.

	longaxis=np.vstack((a,a+vv[0,:]*10))
	import matplotlib.pyplot as plt
	import mpl_toolkits.mplot3d as m3d
	if 0:
		ax = m3d.Axes3D(plt.figure())
		ax.scatter3D(*r.T)
		ax.plot3D(*linepts.T)
		plt.show()

	return longaxis

def findlongvector(r):
	a=r.mean(axis=0)
	uu, dd, vv = np.linalg.svd(r - a)
	# Now vv[0] contains the first principal component, i.e. the direction
	# vector of the 'best fit' line in the least squares sense.

	# Now generate some points along this best fit line, for plotting.

	# I use -7, 7 since the spread of the data is roughly 14
	# and we want it to have mean 0 (like the points we did
	# the svd on). Also, it's a straight line, so we only need 2 points.
	
	return np.vstack((-vv[0,:],vv[0,:]))

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


