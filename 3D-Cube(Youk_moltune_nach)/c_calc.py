import matplotlib.pyplot as plt
import numpy as np
import math

import classes
import random as rd

#Variablen fuer Cgradient#

gamma_ = 7.0                  #degradations constant#
diff_const = 7.0              #diffusions constant#
bruch = float(diff_const/gamma_)
lambda_ = float(math.sqrt(bruch))   #radius of signalcloud of cell#

#Parametersettings#

x=10                              #sets gridsize
n=1.0                             #set chance of cells n probability
place=0.5                         #set on_state cells n probability
C_on=13.0                         #signalconcentration of activ cel#
K =  18.0                         #threshold c#
feedback = 1                      #positiv(1) or negative(0) feedback#
min_cell=5                        #set minimum of cellneigbors for new cellcreation#


ce=0                              #for pos_dic entrykey
pos={}                            #dic for grid and cells in a list [0]=row [1]=col [2] if cell there or not


for row in range(x):
    for col in range(x):

        put=rd.random()
        if put<= n:
            pos[ce]=[row,col,True]     #dict for cell positions#
            ce+=1
        else:
            pos[ce]=[row,col,False]
            ce+=1



r = np.zeros(len(pos))        #list of cellposition in space#
state=np.arange(len(pos))     #list of default state of each cell#
C_i = np.zeros(len(pos))      # concentration of each cell at time i #
C_print=np.zeros(len(pos))    #actual concentrations at position cells + neighbor#
resident=np.zeros(len(pos))   #spot has cell or not#



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
        r[i] = math.sqrt((pos[i][0] - pos[c_num][0]) ** 2 + (pos[i][1] - pos[c_num][1]) ** 2)
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


def calc_expterm(i):
    f_n = 0
    calcDistances(i)
    for j in range(len(r)):          #calc the signal strengh f_n for cell i #
        if j != i:
            f_n = np.exp(-r[j]/lambda_)+f_n
    return f_n


def set_state():                    #function to  set concentration#

    if feedback==1:
        for j, val in enumerate(state):
            if pos[j][2]:
                if val == 1:                #anyone with 1 is activ
                    C_i[j]=C_on             #and concentration is set on c_on for active cell
                    C_print[j]=C_on
                else:
                    C_i[j]=1                #or inactive than set to 1
                    C_print[j]=1


    elif feedback==0:
        for j, val in enumerate(state):
            if pos[j][2]:
                if val == 0:  # anyone with 1 is activ
                    C_i[j] = C_on  # and concentration is set on c_on for active cell
                else:
                    C_i[j] = 1  # or inactive than set to 1

def create(i):                          #function who looks at dircet neighbors to decide if cell comes to life

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

def switch (ci):             # determine cell cis status for next step.#

    on=False                            #boolean to determine cells next state

    #f_n = calc_expterm(ci)
    #print(f_n)
    if feedback==1:           #positiv feedback


        if  C_i[ci]-K >= 0:                         # active state of autonomous cell#
            on=True
            if not pos[ci][2]:
                if create(ci):
                    pos[ci][2]=True
                    on=False


        elif C_i[ci] - K <= 0:                      # deactive state of autonomous cell#
            #print('auto de')
            on = False


        if pos[ci][2]:
            if on:
                state[ci] = 1
                C_i[ci] = C_on
            else:
                state[ci] = 0
                C_i[ci] = 1
        else:
            state[ci] = 0
            C_i[ci] = 0


    elif feedback == 0:  # negative feedback#


        if C_i[ci] - K >= 0:                # deactive state of autonomous cell#

            on = False


        elif C_i[ci] - K <= 0:              # active state of autonomous cell#
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

    for i in pos.keys():

        temp_x.append(pos[i][0])
        temp_y.append(pos[i][1])

    zplot = plt.subplot(3, 2, z)
    zplot.axis([-1, x + 1, -1, x + 1])  # plot of state/cell position

    if not z%2==0:
        if z==1:
            zplot.set_title('Initial concentration')

        elif z==3:
            zplot.set_title('Halftime concentration')

        elif z==5:
            zplot.set_title('End concentration')

        temp_C=np.zeros(len(pos))
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

        for i in pos.keys():
            if pos[i][2]:
                if state[i] == 1:  # if on red dot
                    zplot.plot(pos[i][0], pos[i][1], 'ro')
                else:  # else(if off) blue dot
                    zplot.plot(pos[i][0], pos[i][1], 'bo')

def update(end):

    plt.style.use('ggplot')
    start=0
    stop = end
    time = np.linspace(start, stop, end*10)
    timer=0
    i = 0


    while i < len(pos):             # produce random state

        hp = rd.random()
        if pos[i]:
            if hp <= place:
                state[i] = 1
                i += 1
            else:
                state[i] = 0
                i += 1
        else:
            state[i]= 0
            C_i[i]= 0
            i +=1


    set_state()

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




    print(C_t)
    print(state_t)


print(parameter)
update(1)
plt.show()
