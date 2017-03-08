import numpy as np
"""""
Class that creats pic of concentration gradient of the whole modleroom
"""""
class Gridpic:

    def __init__(self, x, K, c_ary, ddd, pos, C_true):
        # Parametersettings#

        self.pos = pos
        self.x = x
        self.ddd = ddd
        self.K = K
        self.c_ary = c_ary
        self.C_true = C_true
        self.gamma_ = 7.0
        self.diff_const = 7.0
        bruch = float(self.diff_const / self.gamma_)
        self.lambda_ = float(np.sqrt(bruch))


        self.r = np.zeros(len(pos))  # list of cellposition in space#
        self.C_print = np.zeros(len(pos))

    def calcDistances(self,c_num):

        for i in range(len(self.c_ary)):
            self.r[i] = np.sqrt(np.square(self.c_ary[i].xcor - self.pos[c_num][0])
                           + np.square(self.c_ary[i].ycor - self.pos[c_num][1])
                           + np.square(self.c_ary[i].zcor - self.pos[c_num][2]))

        return self.r

    def calc_cval(self):

        for i in range(len(self.pos)):
            c_neighbor = 0.0
            self.calcDistances(i)

            for j in range(len(self.C_true)):  # nachbarzellen_c aufaddieren

                if self.r[j] <= (self.c_ary[j].radius):          #if r[j]< radius they touch so conz should =surface
                    c_neighbor = self.C_true[j] + c_neighbor
                else:
                    c_neighbor = self.C_true[j] \
                                * np.exp(-(self.r[j]) / self.lambda_)\
                                + c_neighbor


            c_value = c_neighbor


            self.C_print[i] = c_value

        return self.C_print
