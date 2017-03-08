import numpy as np


class Cell:
    """
    @type name: str(which bacteria)
    @type radius: float
    @type mid: int
    @type status: bool
    """

    def __init__(self, mid, name, radius, status, xcor, ycor, zcor):
        self.mid = mid
        self.name = name
        self.radius = radius
        #self.molcomp = molcomp
        self.status = status
        self.xcor = xcor
        self.ycor = ycor
        self.zcor = zcor

    @property
    def mid(self):
        return self.__mid

    @mid.setter
    def mid(self, value):
        if not (isinstance(value, float) or isinstance(value, int)):
            raise Exception("radius must be numeric")
        else:
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

    @property
    def x_cor(self):
        return self.__x_cor

    @x_cor.setter
    def x_cor(self, value):
        if not (isinstance(value, float) or isinstance(value, int)):
            raise Exception("coordinate must be numeric")
        else:
            self.__x_cor = value

    @property
    def y_cor(self):
        return self.__y_cor

    @y_cor.setter
    def y_cor(self, value):
        if not (isinstance(value, float) or isinstance(value, int)):
            raise Exception("coordinate must be numeric")
        else:
            self.__y_cor = value

    @property
    def z_cor(self):
        return self.__z_cor

    @z_cor.setter
    def z_cor(self, value):
        if not (isinstance(value, float) or isinstance(value, int)):
            raise Exception("coordinate must be numeric")
        else:
            self.__z_cor = value

    def __repr__(self):
        return self.name

class C_grad:

    def __init__(self, c_num,x,n,place,C_on,K,feedback,c_ary,ddd):
        #Parametersettings#
        self.c_num = c_num
        self.x = x
        self.n = n
        self.ddd=ddd
        self.place = place
        self.C_on = C_on
        self.K = K
        self.feedback = feedback
        self.c_ary = c_ary
        self.gamma_ = 4.0
        self.diff_const = 7.0
        bruch = float(self.diff_const / self.gamma_)
        self.lambda_ = float(np.sqrt(bruch))
        self.P=2000                         #particlenumber of production

        self.r = np.zeros(c_num)            #list of cellposition in space#
        self.state = np.zeros(c_num)  # list of default state of each cell#
        self.C_i = np.zeros(c_num)          #concentration of each cell at time i #
        self.C_print = np.zeros(c_num)      #actual concentrations at position cells + neighbor#


    @property
    def c_num(self):
        return self.__c_num
    @c_num.setter
    def c_num(self, value):
        if not (isinstance(value, float) or isinstance(value, int)):
            raise Exception("There must be cells")
        else:
            self.__c_num = value

    def calcDistances(self,c_num):

        for i in range(len(self.r)):
            self.r[i] = np.sqrt(np.square(self.c_ary[i].xcor - self.c_ary[c_num].xcor)
                           + np.square(self.c_ary[i].ycor - self.c_ary[c_num].ycor)
                           + np.square(self.c_ary[i].zcor - self.c_ary[c_num].zcor))

        return self.r

    def calc_cval(self,step, C_t):

        for i in range(len(self.C_i)):
            c_neighbor = 0.0
            self.calcDistances(i)

            for j in range(len(self.C_i)):  # nachbarzellen_c aufaddieren

                conz_j=self.catcher_c(step,C_t,j)

                if j == i:
                    continue

                if self.c_ary[j].name == 'lazy':           #produce nothing so no calculation
                    continue

                if j != i:

                    if self.r[j] <= self.c_ary[i].radius:          #if r[j]< radius they touch so conz should =surface
                        c_neighbor += conz_j
                    else:
                        c_neighbor += conz_j \
                                      * (self.c_ary[i].radius / self.r[j]) \
                                      * np.exp(-(self.r[j] - self.c_ary[i].radius) / self.lambda_)

            conz_i = self.catcher_c(step, C_t, i)
            c_value = c_neighbor + conz_i

            if self.c_ary[i].name == 'lazy':     #have their default c
                c_value = c_neighbor

            self.C_i[i] = c_value
            self.C_print[i] = c_value

        return self.C_i

    def catcher_c(self,step,C_t,j):
        try:
            conz = C_t[step][j]

        except IndexError:
            if self.c_ary[j].status == False:
                if self.feedback == 1:
                    conz = 1                    #inactive than set to 1
                elif self.feedback == 0:        #active set to C_on
                    conz = self.C_on
            else:
                if self.feedback == 1:
                    conz = self.C_on
                elif self.feedback == 0:
                    conz = 1
        if self.ddd:
            conz=self.conz_r(j)

        return conz

    def conz_r(self,i):                   #calc c for sphere cell

        top=(self.P*self.c_ary[i].radius)
        down1=(3*np.sqrt(self.diff_const*self.gamma_)*self.c_ary[i].radius)
        down2=self.diff_const*3
        c_r=top/(down1+down2)

        return c_r

    def ini_cell_c(self):

        for i in range(len(self.c_ary)):  # produce random state

            if self.c_ary[i].status:
                self.state[i] = 1
                self.set_c(i)

            else:
                self.state[i] = 0
                self.set_c(i)

        return (self.C_i , self.state)

    def set_c(self,p):

        if self.c_ary[p].status == False:

            if self.feedback == 1:
                self.C_i[p] = 1  # or inactive than set to 1
                self.C_print[p] = 1

            elif self.feedback == 0:
                self.C_i[p] = self.C_on  # and concentration is set on c_on for active cell
                self.C_print[p] = self.C_on

        elif self.c_ary[p].status == True:

            if self.feedback == 1:
                self.C_i[p] = self.C_on  # and concentration is set on c_on for active cell
                self.C_print[p] = self.C_on

            elif self.feedback == 0:
                self.C_i[p] = 1  # or inactive than set to 1
                self.C_print[p] = 1

    def switch(self):  # determine cell cis status for next step.#

        for ci in range(len(self.C_i)):

            on = False  # boolean to determine cells next state

            if self.ddd:

                top1 = (1/self.C_on)
                top2 = (self.P * self.c_ary[ci].radius)
                down1 = (3 * np.sqrt(self.diff_const * self.gamma_) * self.c_ary[ci].radius)
                down2 = self.diff_const * 3
                fac=top1*top2/(down1+down2)

                K_r=self.K*fac
            else:
                K_r = self.K

            if self.feedback == 1:  # positiv feedback

                if self.C_i[ci] - K_r > 0:  # active state of autonomous cell#
                    on = True

                elif self.C_i[ci] - K_r <= 0:  # deactive state of autonomous cell#
                    on = False

                #if off stays off if on cell can change
                if on:
                    self.state[ci] = 1
                    self.C_i[ci] = self.C_on
                else:
                    self.state[ci] = 0
                    self.C_i[ci] = 1

            elif self.feedback == 0:  # negative feedback#

                if self.C_i[ci] - K_r > 0:  # deactive state of autonomous cell#
                    on = False

                elif self.C_i[ci] - K_r <= 0:  # active state of autonomous cell#
                    # print('auto ac')
                    on = True

                if on:
                    self.state[ci] = 1
                    self.C_i[ci] = self.C_on
                else:
                    self.state[ci] = 0

                    self.C_i[ci] = 1

        return self.C_i, self.state

class Molecule:
    """

    @type name: str
    @type mass: float
    """

    def __init__(self, mid, name, mass=0, count=0):
        self.mid = mid
        self.name = name
        self.mass = mass
        self.count = count

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

    @property
    def count(self):
        return self.__count

    @count.setter
    def count(self, value):
        if not (isinstance(value, float) or isinstance(value, int)):
            raise Exception("mass must be numeric")
        else:
            self.__count = value

    def __repr__(self): #string "self.name"		#print(list(object))
        return self.name

class BioMolecule:
    """
    A generic molecule that has basic attributes like name and
    mass.
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

    def __repr__(self): #string "self.name"
        return self.name

class BioMoleculeCount:        #new variable: number of molecules of this type


    def __init__(self, mid, name, count=0):
        super(BioMolecule).__init__(mid, name)
        self.count = count

    @property
    def count(self):
        return self.__count

    @count.setter
    def count(self, value):
        self.__count = value