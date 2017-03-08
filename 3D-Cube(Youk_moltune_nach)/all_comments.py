import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib
import random as rd
import hl
import shaft
#Variablen fuer Cgradient#

gamma_ = 1                   #degradations constant#
diff_const = 1             #diffusions constant#
bruch = diff_const/gamma_
lambda_ = math.sqrt(bruch)   #radius of signalcloud of cell#
#Parametersettings

x=10                         #sets gridsize
n=1                         #sets chance of cells (1/n)
pos={}
ce=0

for a in range(x):
    for b in range(x):
        put=rd.randint(1, n)
        if put== n:
            pos[ce]=[a,b] # dict for cell positions#
            ce+=1

feedback = 1                 #positiv(1) or negative(0) feedback#

r = np.zeros(len(pos))            #list of cellposition in space#
state=np.arange(len(pos))          #list of default state of each cell#
C_i = np.zeros(len(pos))          # concentration of each cell at time i #
C_print=np.zeros(len(pos))

C_on=1.5                      #signalconcentration of activ cel#
K =  7                     #threshold c#

#outputdata

C_t=[]                       # concentrations at time step#
state_t=[]        #lists all individual cellstates at timestep#

parameter = 'C_on = ' + str(C_on) \
            + ' K = ' + str(K) \
            + ' n = ' + str(n) \
            + ' gamma  = ' + str(gamma_) \
            + ' diffconst =' + str(diff_const) \
            + ' feedback = ' + str(feedback) \
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
    #print('c_value :'+ str(c_neighbor))
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
            if val == 1:                #anyone with 1 is activ
                C_i[j]=C_on             #and concentration is set on c_on for active cell
                C_print[j]=C_on
            else:
                C_i[j]=1                #or inactive than set to 1
                C_print[j]=1

    elif feedback==0:
        for j, val in enumerate(state):
            if val == 0:  # anyone with 1 is activ
                C_i[j] = C_on  # and concentration is set on c_on for active cell
            else:
                C_i[j] = 1  # or inactive than set to 1




def neighbor():                    #function to sum the current activated neighboring cells #
    acti = 0.0

    for j, val in enumerate(state):
        if val == 1:      #anyone else whos activ is counted
            acti += 1

    return acti


def switch (step,ci, acti):             # determine cell cis status for next step.#

    on=False                            #boolean to determine cells next state
    #quorum=False
    #all=float(len(state))

    #if (acti/all) >= 0.5:            #quorum sensing possibility if more active than cell also active#
        #print(acti)
        #quorum=True


    #f_n = calc_expterm(ci)
    #print(f_n)


    #print('Soviele active drum rum '+str(ci)+' '+str(acti))
    #print('f_n ='+str(f_n))

    #A_0 = 1 + f_n - K
    #A_n = C_on + (1-K) / f_n
    #+ (1 + (len(state) - acti) * f_n) / f_n  # activation  threshold

    #D_0 = C_i[ci] - K + f_n
    #D_n = C_on - K / (f_n + 1)
    # + (1 + (acti) * f_n) / (1 + f_n)  # deactivation threshold


    if feedback==1:           #positiv feedback

        #print('A_n ='+str(A_n))
        #print('A_0 =' + str(A_0))
        #print('D_0 ='+str(D_0))
        #print('D_n =' + str(D_n))
        #print('C_i =' + str(C_i[ci]))
        #print(str(C_i[ci]-K))

        if  C_i[ci]-K >= 0:                         # active state of autonomous cell#

            # if inbetween bistable-state of autonomous cell#

            #if state[ci] == 0:
                #on=False

           # else:
            on=True

            #on=True


        elif C_i[ci] - K <= 0:                      # deactive state of autonomous cell#
            #print('auto de')
            on = False


        #if C_i[ci] >= D_0 and C_i[ci] <= A_n:  # if inbetween bistable-state of autonomous cell#
            #print('auto bi')
            #if state[ci] == 0:
                #on=False

            #else:
                #on=True


        #if  C_i[ci]<= A_n and C_i[ci]>= A_0 :   #if true than active in next state cause of neighbor
            #print('ne ac')
            #on=True

        #if  C_i[ci]>= D_n and C_i[ci]<=D_0:   # if true than deactive in next state because of neighbor
            #print('ne de')
            #on=False

        #if C_i[ci] <= D_0 and C_i[ci] >= A_n:  # if in between quorum cell#
            #print('ne qu')
            #if quorum:
                #print('ne qu')
                #on = True

            #else:
                #on = False

        if on:
            state[ci] = 1
            C_i[ci] = C_on
        else:
            state[ci] = 0
            C_i[ci] = 1




    elif feedback == 0:  # negative feedback#

        #print('A_n ='+str(A_n))
        #print('D_0 ='+str(D_0))
        #print('C_i =' + str(C_i[ci]))

        if C_i[ci] - K >= 0:                # deactive state of autonomous cell#
            #print('auto de')
            on = False


        elif C_i[ci] - K <= 0:              # active state of autonomous cell#
            #print('auto ac')
            on = True


        #elif C_i[ci] >= D_0 and C_i[ci] <= A_n:  # if inbetween flipflop of autonomous cell#
            #print('auto flip')
            #if state[ci] == 0:
                #on = True

            #else:
                #on = False

        #elif C_i[ci] >= A_n and C_i[ci] <= A_0:  # if true than deactive in next state cause of neighbor
            #print('ne de')
            #on = False

        #elif C_i[ci] >= D_n and C_i[ci] <= D_0:  # if true than active in next state
            #print('ne ac')
            #on = True

        #elif C_i[ci] <= D_0 and C_i[ci] >= A_n:  # if in between quorum cell#
            #print('ne qu')

            #if quorum:

                #on = False

            #else:
                #on = True

        if on:
            state[ci] = 1
            C_i[ci] = C_on
        else:
            state[ci] = 0
            C_i[ci] = 1

def figureprint(z):
    print(C_print)
    temp_x = []
    temp_y = []

    for i in pos.keys():
        temp_x.append(pos[i][0])
        temp_y.append(pos[i][1])

    z = plt.subplot(3, 1, z)
    z.axis([-1, x + 1, -1, x + 1])  # plot of state/cell position

    if z==1:
        z.set_title('Inital state')

    elif z==2:
        z.set_title('Halftime state')

    elif z==3:
        z.set_title('End state')

    temp_C=np.zeros(len(pos))
    for i in range(len(C_print)):
        temp_C[i]=C_print[i]/max(C_print)


    z.scatter(temp_x, temp_y, c=temp_C, cmap='bwr')


    for i in pos.keys():
        if state[i] == 1:  # if on red dot
            plt.plot(pos[i][0], pos[i][1], 'ro')
        else:  # else(if off) blue dot
            plt.plot(pos[i][0], pos[i][1], 'bo')


def update(end):

    plt.style.use('ggplot')
    start=0
    stop = end
    time = np.linspace(start, stop, end*10)
    timer=0
    i = 0
    #50_50 positioning
    while i < len(pos):
        if i <= len(pos) / 2:
            state[i] = 1  # rd.randint(0, 1)  # produce random state
            i += 1
        if i >= len(pos) / 2:
            state[i] = 0  # rd.randint(0, 1)  # produce random state
            i += 1

    while i < len(pos):
        state[i] = rd.randint(0, 1)  # produce random state
        i += 1
    #stat_p=str(state)


    for step in time:

        set_state()

        if timer==0:                             #saves initial status of system
            C_t.append(list(C_i))
            state_t.append(list(state))


            figureprint(1)



        for i in range(len(C_i)):                 #calc actual c for cell
            calc_cval(timer,i)

        acti = neighbor()

        for j in range(len(C_i)):                 #determine cells behaviour
            switch(timer,j,acti)


        state_t.append(list(state))              #saves status
        C_t.append(list(C_i))                    #saves set c of cells

        if timer * 2 == len(time):
            figureprint(2)

        if timer == len(time) - 2:
            figureprint(3)
            fig.savefig(parameter, dpi=600, format='pdf', bbox_inches='tight')
        timer += 1



    print(C_t)
    print(state_t)


print(parameter)
update(1)

plt.show()

# Dimensions
nx, ny, nz = x, x, 1
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
dis=np.zeros((nx,ny))

def calcDis_vis(x,y):

    for i in range(nx):
        for j in range(ny):
            dis[i][j] = math.sqrt((i - x) ** 2 + (j - y) ** 2)


    return dis


def distribution_fct(i):

    c_neighbor = 0
    calcDistances(i)
    for j in range(len(C_i)):  # nachbarzellen_c aufaddieren
        if j != i:
            c_neighbor = C_i[j] * np.exp(-r[j]/lambda_) + c_neighbor

    c_value = c_neighbor + C_i[i]
    # print('c_value :'+ str(c_neighbor))

    return c_value
#p=0
#for i in range(0,nx):
    #for j in range(0,ny):
        #distribution[i][j][0] = distribution_fct(p)
        #p+=1
        #print "x:" +  str(x[i])
        #print "y:" + str(y[j])
        #print "value:" + str(distribution[i][j][0])



