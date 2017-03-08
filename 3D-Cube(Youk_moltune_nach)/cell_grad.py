import matplotlib.pyplot as plt
import numpy as np
import random as rd
import gridlib as gl

#Variablen fuer Cgradient#

gamma_ = 7.0                            #degradations constant#
diff_const = 1.0                        #diffusions constant#
bruch = float(diff_const/gamma_)
lambda_ = float(np.sqrt(bruch))       #radius of signalcloud of cell#

#Parametersettings#

x = gl.x                                  #sets gridsize
n = gl.n                                  #set chance of cells n probability
place = gl.place                          #set on_state cells n probability
C_on = gl.C_on                            #signalconcentration of activ cell#
K =  gl.K                               #threshold c#
feedback = gl.feedback                  #positiv(1) or negative(0) feedback#
#min_cell=5                             #set minimum of cellneigbors for new cellcreation#
c_ary=gl.c_ary

r = np.zeros(len(gl.c_ary))             #list of cellposition in space#
state = np.zeros(len(gl.c_ary))         #list of default state of each cell#
C_i = np.zeros(len(gl.c_ary))           #concentration of each cell at time i #
C_print = np.zeros(len(gl.c_ary))       #actual concentrations at position cells + neighbor#

#outputdata#

C_t=[]                        # concentrations at time step#
state_t=[]                    #lists all individual cellstates at timestep#

parameter = 'C_on = ' + str(C_on) \
            + ' K = ' + str(K) \
            + ' n = ' + str(n) \
            + ' gamma  = ' + str(gamma_) \
            + ' diffconst =' + str(diff_const) \
            + ' feedback = ' + str(feedback) \
            + ' on_per_c = '+ str(place)\
            + ' x = ' + str(x)\
            +'.pdf'

fig = plt.figure()
fig.suptitle('Cell position and state', fontsize=16)



def calcDistances(c_num):
    for i in range(len(r)):
        r[i] = np.sqrt(  np.square(c_ary[i].xcor - c_ary[c_num].xcor)
                       + np.square(c_ary[i].ycor - c_ary[c_num].ycor)
                       + np.square(c_ary[i].zcor - c_ary[c_num].zcor) )
        #print(r[i])

    return r



def calc_cval(step,i):
    c_neighbor = 0.0
    calcDistances(i)
    #print(r)
    for j in range(len(C_i)):                #nachbarzellen_c aufaddieren
        if j!=i:
            c_neighbor = C_t[step][j] * np.exp(-r[j] / lambda_) + c_neighbor

    c_value = c_neighbor + C_t[step][i]

    C_i[i]=c_value
    C_print[i]=c_value



"""""
def create(i):
#function who looks at dircet neighbors to decide if cell comes to life

    counter=0
    life=False

    if i-1>0:
        if pos[i-1] and state[i-1]:
            counter+=1
    if i+1<x*x:
        if pos[i+1] and state[i-1]:
            counter+=1
    if i-x >= 0:
        if pos[i-x]and state[i-x] :
            counter+=1
    if i+x < x*x:
        if pos[i+x]and state[i+x]:
            counter+=1

    if counter>=min_cell:
        life=True

    return life
"""""
def switch (ci):                        #determine cell cis status for next step.#

    on=False                            #boolean to determine cells next state

    if feedback==1:                     #positiv feedback

        if  C_i[ci]-K >= 0:             #active state of autonomous cell#
            on=True
            #if not pos[ci][2]:
                #if create(ci):
                    #pos[ci][2]=True
                    #on=False
        elif C_i[ci] - K <= 0:          #deactive state of autonomous cell#

            on = False

        if c_ary[ci].status:            #if off stays off if on cell can change
            if on:
                state[ci] = 1
                C_i[ci] = C_on
            else:
                state[ci] = 0
                C_i[ci] = 1
        else:
            state[ci] = 0
            C_i[ci] = 0


    elif feedback == 0:                     #negative feedback#


        if C_i[ci] - K >= 0:                #deactive state of autonomous cell#
            on = False

        elif C_i[ci] - K <= 0:              #active state of autonomous cell#
            #print('auto ac')
            on = True

        if on:
            state[ci] = 1
            C_i[ci] = C_on
        else:
            state[ci] = 0
            C_i[ci] = 1


def figureprint (z):

    #print(C_print)
    temp_x = []
    temp_y = []

    for i in c_ary.keys():

        temp_x.append(c_ary[i].xcor)
        temp_y.append(c_ary[i].ycor)

    zplot = plt.subplot(3, 2, z)
    zplot.axis([-1, x + 1, -1, x + 1])  # plot of state/cell position

    if not z%2==0:
        if z==1:
            zplot.set_title('Initial concentration')

        elif z==3:
            zplot.set_title('Halftime concentration')

        elif z==5:
            zplot.set_title('End concentration')

        temp_C=np.zeros(len(c_ary))
        for i in range(len(C_print)):
            temp_C[i]=C_print[i]/K


        zplot.scatter(temp_x, temp_y, c=temp_C, cmap='bwr')
    if  z % 2 == 0:
        if z == 2:
            zplot.set_title('Initial state')

        elif z == 4:
            zplot.set_title('Halftime state')

        elif z == 6:
            zplot.set_title('End state')

        for i in c_ary.keys():
            if c_ary[i].status:
                if state[i] == 1:  # if on red dot
                    zplot.plot(c_ary[i].xcor, c_ary[i].ycor, 'ro')
                else:  # else(if off) blue dot
                    zplot.plot(c_ary[i].xcor, c_ary[i].ycor, 'bo')

def ini_cell_c():

    for i in range(len(c_ary) + 1):           #produce random state

        hp = rd.random()
        if c_ary[i].status:
            if hp <= place:
                state[i] = 1
                set_c(i)
            else:
                state[i] = 0
                set_c(i)
        else:
            state[i] = 0
            set_c(i)


def set_c(p):

     if state(p) == 0:
         if feedback == 1:
             C_i[p] = 1                     #or inactive than set to 1
             C_print[p] = 1
         elif feedback == 0:
             C_i[p] = C_on                  #and concentration is set on c_on for active cell
             C_print[p] = C_on

     elif state(p) == 1:
         if feedback == 1:
             C_i[p] = C_on  # and concentration is set on c_on for active cell
             C_print[p] = C_on
         elif feedback == 0:
             C_i[p] = 1  # or inactive than set to 1
             C_print[p] = 1

def update(end):

    plt.style.use('ggplot')
    start=0
    stop = end
    time = np.linspace(start, stop, end*10)
    timer=0
    i = 0



    ini_cell_c()

    for step in time:

        if timer==0:                             #saves initial status of system
            C_t.append(list(C_i))
            state_t.append(list(state))
            figureprint(1)
            figureprint(2)

        for i in range(len(C_i)):                 #calc actual c for cell
            calc_cval(timer,i)


        for j in range(len(C_i)):                 #determine cells behaviour
            switch(j)


        state_t.append(list(state))              #saves status

        C_t.append(list(C_i))                    #saves set c of cells

        if timer * 2 == len(time):
            figureprint(3)
            figureprint(4)

        if timer == len(time) -2:
            figureprint(5)
            figureprint(6)

            fig.savefig(parameter, dpi=600, format='pdf', bbox_inches='tight')
        timer += 1

print(parameter)
update(1)
plt.show()
