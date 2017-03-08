import numpy as np

class pm:
    # Parametersettings#
    x = 10  # sets gridsize
    n = 1.0  # set chance of cells n probability
    place = 0.6  # set on_state cells n probability
    C_on = 19.0  # signalconcentration of activ cel#
    K = 18.0  # threshold c#
    feedback = 1  # positiv(1) or negative(0) feedback#

    radius = 1  # set intial radius of cell
    maxrad = 2 * radius
    maxvolume = 2 * 4 * np.pi * radius ** 3  # max r which cell can have before division#

    k = 1.0  # factor how close neighboring cells can be 1.0 just touching
    g_rate = 0.2  # growthrate of radius

    F_g = 0.002  # gravity constant
    threeD = 1  # 3d or 2d

    thres_skr = 0.3  # threshold for skrinking

    cell_age = 0
    maxcell_age = 6

    # Dimensions
    nx, ny, nz = x, x, 1
    if threeD:
        nz = x
    lx, ly, lz = 10.0, 10.0, 0.1
    if threeD:
        lz = 10.0
    dx, dy, dz = lx / nx, ly / ny, lz / nz

    ncells = nx * ny * nz
    npoints = (nx + 1) * (ny + 1) * (nz + 1)

    # Coordinates
    x_c = np.arange(0, lx + dx, dx, dtype='float64')
    y_c = np.arange(0, ly + dy, dy, dtype='float64')
    z_c = np.arange(0, lz + dz, dz, dtype='float64')

print pm.maxvolume