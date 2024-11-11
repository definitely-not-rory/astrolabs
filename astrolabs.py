from imports import *
from data_handling import get_data, raw_plot
from fourier import fourier_fitting
from sinusoidal import sin_fitting
from plotgrids import plot_grid
from alias import alias
from colour_plot import colour_plot

df = pd.read_csv('mcmaster.txt')

arrays=df.to_numpy()

objs=next(os.walk('.'))[1]

objs=objs[2:-1]
display_grid=input('Show Grid (Y/N): ')

if display_grid=='Y' or display_grid=='y':
    plot_grid(objs)

obj=input('Select Cepheid: ')

period=0

for i in arrays:
    if i[0]==obj:
        period=i[1]
    
print('~~~ '+obj+' ~~~\nLiterature Period Value: '+str(period))
try:
    fourier_fitting(obj,period,2,5,True)
except:
    print('Not Enough Data To Use Fourier Modelling')



'''
RORY TO DO:

change period % bound

grid plot of fourier and folded curves

residuals

'''