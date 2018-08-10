import numpy as np
import scipy.integrate as integrate

# Integration of variable bins using
# Rectacgles (rect), Trapezoids (trap)
def int_vbin(y,x,xlow,xhigh,int_type='trap'):
    val = 0.
    if (int_type == 'rect'):
        for i,h in enumerate(y):
            val = val + (xhigh[i]-xlow[i])*h
    elif (int_type == 'trap'): # same result as np.trapz(y,x)
        for i in range(len(y)-1):
            val = val + (x[i+1]-x[i])*(y[i+1]+y[i])/2.
    return val
