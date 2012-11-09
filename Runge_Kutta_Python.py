from math import *
from numpy import *
from matplotlib import *
from pylab import *

def eulermethod(the_deriv,y0,t0,dt,n):
    y_list = [y0]
    t_list = [t0]
    for i in range(0,n):
        m = the_deriv(t0,y0)
        y0 = y0 + m*dt
        t0 = t0 + dt
        y_list.append(y0)
        t_list.append(t0)
        y0 = y0
        t0 = t0
    return([t_list,y_list])



#Uses fourth order Runge-Kutta algorithm to calculate the trajectory of a derivative given the derivative as a python function, initial values y0 and t0, time step dt, and the number of steps n. Returns vector of the t values vector and the y values vector
def rk4(the_deriv,t0,y0,dt,n):
    y_list = [y0]
    t_list = [t0]
    #iterate through the rk4 calculation for n number of steps
    for i in range(0,n):
        #calculate the k-values for rk4
        k1 = dt*the_deriv(t0,y0)
        k2 = dt*the_deriv((t0+dt/2),(y0+k1/2))
        k3 = dt*the_deriv((t0+dt/2),(y0+k2/2))
        k4 = dt*the_deriv((t0+dt),(y0+k3))
        #calculate y1 based on k-values, add dt to t, append to the list
        y1 = y0 +k1/6+k2/3+k3/3+k4/6
        t1 = t0+dt
        y_list.append(y1)
        t_list.append(t1)

        #reset the values for the next cycle
        y0 = y1
        t0 = t1
    return([t_list,y_list])




###Compare results with these tests
##the derivative function
def the_deriv(t,y):
    out = -2.3*y
    return(out)
#
##initial values be sure to use floats (i.e. 1=1.0)
y0=1.0
t0=0.0
tmax = 15
dt=.0125
n=int(tmax/dt)


###Euler Metohd Test
#test = eulermethod(the_deriv,y0,t0,dt,n)
#t = test[0]
#y = test[1]
#plot(t,y)
#
#RK4 Test
test2 = rk4(the_deriv,t0,y0,dt,n)
t = test2[0]
y = test2[1]
print(t)
print(y)
#plot(t,y)
#
###Actaul answer Test
#t = arange(0,tmax,dt)
#y = exp(-2.3*t)
#plot(t,y)
#
#plt.legend(['Euler', 'RK4', 'Actual'], loc='upper right')
#
#show()




