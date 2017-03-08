#! /usr/bin/env python

import hl
import numpy as np

# Dimensions
nx, ny, nz = 40, 40, 1
lx, ly, lz = 10.0, 10.0, 0.1
dx, dy, dz = lx/nx, ly/ny, lz/nz

ncells = nx * ny * nz
npoints = (nx + 1) * (ny + 1) * (nz + 1)

# Coordinates
x = np.arange(0, lx + dx, dx, dtype='float64')
y = np.arange(0, ly + dy, dy, dtype='float64')
z = np.arange(0, lz + dz, dz, dtype='float64')

# Variables
point_data = np.zeros((nx + 1, ny + 1, nz + 1))
#cell_data = np.zeros((nx, ny, nz))
distribution = np.zeros((nx, ny, nz))



D = 10.0
k = 0.5
lam = np.sqrt(D/k)
cell_radius = 2.0



def distribution_fct(x,y):
  x0 = 5.0
  y0 = 5.0
  r = max(cell_radius,np.sqrt((x-x0)**2 + (y-y0)**2))
  value = np.exp(-r/lam)/r
  return value

for i in range(0,nx):
    for j in range(0,ny):
      distribution[i][j][0] = distribution_fct(x[i],x[j])
      #print "x:" +  str(x[i])
      #print "y:" + str(y[j])
      #print "value:" + str(distribution[i][j][0])	



hl.gridToVTK("./test", x, y, z, cellData = {"distribution" : distribution}, pointData = {"temp" : point_data})