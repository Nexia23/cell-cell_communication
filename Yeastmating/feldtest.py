import numpy as np
threeD=1
x=10
# Dimensions
nx, ny, nz = x, x, 1
if threeD:
    nz = x
lx, ly, lz = 10.0, 10.0, 0.1
if threeD:
    lz = 10.0
dx, dy, dz = lx/nx, ly/ny, lz/nz

ncells = nx * ny * nz
npoints = (nx + 1) * (ny + 1) * (nz + 1)

# Coordinates
x_c = np.arange(0, lx + dx, dx, dtype='float64')
y_c = np.arange(0, ly + dy, dy, dtype='float64')
z_c = np.arange(0, lz + dz, dz, dtype='float64')


print(len(x_c))