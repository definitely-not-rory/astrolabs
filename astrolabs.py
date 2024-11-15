from imports import *
from data_handling import get_data, raw_plot
from fourier import fourier_fitting
from sinusoidal import sin_fitting
from alias import alias
from colour_plot import colour_plot
from fourier_grid import plot_fourier_grid
from fourier_bounds import *

df = pd.read_csv('mcmaster.txt')

arrays=df.to_numpy()

objs=next(os.walk('.'))[1]

objs=objs[2:-1]
periods=[]

for j in objs:
    for i in arrays:
        if i[0]==j:
            periods.append(i[1])
        

obj=input('Select Cepheid: ')

period=0

for i in arrays:
    if i[0]==obj:
        period=i[1]
    
print('~~~ '+obj+' ~~~\nLiterature Period Value: '+str(period))

'''
RORY TO DO:

chi on number of obs

chi on number of fourier modes

fourier pages

change plots designs

plot curves +/- error in period

residuals

'''
plot_fourier_grid(objs,periods,True)
fourier_fitting(obj,period,2,5,True,True,20)
bounds_plot(obj,period)