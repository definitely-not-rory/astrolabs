from imports import *
from data_handling import get_data, raw_plot
from fourier import fourier_fitting
from sinusoidal import sin_fitting
from plotgrids import plot_grid
from alias import alias

df = pd.read_csv('mcmaster.txt')

arrays=df.to_numpy()

obj='ry_cas'

objs=next(os.walk('.'))[1]

objs=objs[2:-1]

period=0

for i in arrays:
    if i[0]==obj:
        period=i[1]
    
print('~~~ '+obj+' ~~~\nLiterature Period Value: '+str(period))

#plot_grid(objs)

alias(obj,period,10,15,300)

'''
RORY TO DO:

jackknifing

colour plot of odd and even coefficients on chi_sq for fourier

plot fitted function on folded data

'''