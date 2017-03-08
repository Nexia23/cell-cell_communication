
import random as rd
import numpy as np
import math

class cell:
    """
    @type name: str(which bacteria)
    @type radius: float
    @type molcomp: dict
    @type status: bool
    """

    def __init__(self, mid, name, radius, status):
        self.mid = mid
        self.name = name
        self.radius = radius
        #self.molcomp = molcomp
        self.status = status


    @property
    def mid(self):
        return self.__mid

    @mid.setter
    def mid(self, value):
        self.__mid = value

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, value):
        self.__name = value

    @property
    def radius(self):
        return self.__radius

    @radius.setter
    def radius(self, value):
        if not (isinstance(value, float) or isinstance(value, int)):
            raise Exception("radius must be numeric")
        else:
            self.__radius = value

    @property
    def status(self):
        return self.__status

    @status.setter
    def status(self,value):
        self.__status = value


    def __repr__(self):
        return self.name

def border(h,b):

    fenc = b + x * h

    if fenc - radius < 0:
        return False
    if fenc + radius > (x) ** 2:
        return False

    if h-radius < 0 :
        return False
    if b-radius< 0 :
        return False

    else:
        return True




def occupy(m):                              #cellplacement

    if isinstance(pos[m][2], cell):         #check if cell there safety
        #print('cell,ho')
        rad = pos[m][2].radius

        for h in range(-rad,rad+1):         # go along the (row)
            for b in range(-rad,rad+1):     #outer quadrate of the circle r=rad (col)

                poi = m + h*x + b  # position in the chain

                eq = math.sqrt((pos[poi][0]-pos[m][0]) ** 2 + (pos[poi][1]-pos[m][1]) ** 2)

                if eq <= rad:               #if smaller or same than position inside the cell

                    if pos[poi][2]== False:
                        pos[poi][2] = str(pos[m][2].mid)

                elif poi==m:
                    pass




def checkneighbor(p):                   # check if can be placed

    ind = True

    if p >= (x+1)**2 or p < 0:
        return False

    if isinstance(pos[p][2], basestring):
        #print('str')
        return False

    if isinstance(pos[p][2], cell):
        #print('cell')
        return False

    rad = radius

    for h in range(-rad, rad + 1):                      # go along the (row)
        for b in range(-rad, rad + 1):                  # outer quadrate of the circle r=rad (col)

            poi = p + h * x + b                         # position in the chain
            if poi >= (x + 1) ** 2 or poi < 0:
                return False

            eq = math.sqrt((pos[poi][0] - pos[p][0]) ** 2 + (pos[poi][1] - pos[p][1]) ** 2)

            if eq <= rad:                               # if smaller or same than position inside the cell

                if isinstance(pos[poi][2], cell):
                    ind = False

                if isinstance(pos[poi][2], basestring):
                    ind = False

    return ind


def move(p):
    pass
    #print('shiiit'+ str(p))

#Parametersettings#

x=5                              #sets gridsize
n=1.0                             #set chance of cells n probability
place=0.5                         #set on_state cells n probability
C_on=13.0                         #signalconcentration of activ cel#
K =  18.0                         #threshold c#
feedback = 0                      #positiv(1) or negative(0) feedback#
min_cell=5                        #set minimum of cellneigbors for new cellcreation#
radius=2

cc=0                              #for pos_dic entrykey
pos={}                            #dic for grid and cells in a list [0]=row [1]=col [2] if cell there or not

# Dimensions
nx, ny, nz = x, x, 1
lx, ly, lz = 10.0, 10.0, 0.1
dx, dy, dz = lx/nx, ly/ny, lz/nz

ncells = nx * ny * nz
npoints = (nx + 1) * (ny + 1) * (nz + 1)

# Coordinates
xs = np.arange(0, lx + dx, dx, dtype='float64')
y = np.arange(0, ly + dy, dy, dtype='float64')
z = np.arange(0, lz + dz, dz, dtype='float64')

r = np.zeros(ncells)        #list of cellposition in space#
state=np.arange(ncells)     #list of default state of each cell#
C_i = np.zeros(ncells)      # concentration of each cell at time i #
C_print=np.zeros(ncells)    #actual concentrations at position cells + neighbor#
resident=np.zeros(ncells)   #spot has cell or not#

def initialize():                       #creats grid on which cells placed
    ce = 0
    for row in range(x+1):
        for col in range(x+1):
            pos[ce] = [row, col, False]
            ce += 1


initialize()

for row in range(x+1):
    for col in range(x+1):
        put=rd.random()

        if border(row,col):

            if checkneighbor(cc):

                if put<= n :

                    pos[cc]= [row,col,cell(cc,'ecoli ', radius, True)]     #dict for cell positions#
                    #print('Hey'+str(cc))
                    occupy(cc)
                    #print(pos)

                    cc+=1
                else:
                    cc+=1

            else:
                cc+=1



        else:
            cc+=1


print(pos)


#print(pos)

class Molecule:
    """

    @type name: str
    @type mass: float
    """

    def __init__(self, mid, name, mass=0):
        self.mid = mid
        self.name = name
        self.mass = mass

    @property
    def mid(self):
        return self.__mid

    @mid.setter
    def mid(self, value):
        self.__mid = value

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, value):
        self.__name = value

    @property
    def mass(self):
        return self.__mass

    @mass.setter
    def mass(self, value):
        if not (isinstance(value, float) or isinstance(value, int)):
            raise Exception("mass must be numeric")
        else:
            self.__mass = value

    def __repr__(self): #string "self.name"		#print(list(object))
        return self.name