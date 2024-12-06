from imports import *
from data_handling import get_data, raw_plot, get_bv
from fourier import fourier_fitting
from sinusoidal import sin_fitting
from alias import alias
from colour_plot import colour_plot
from fourier_grid import plot_fourier_grid
from fourier_bounds import *
from modes_plot import *
from overfitting import *
from chionobs import *
from periodperiod import *
from aperture_photometry_plot import *


df = pd.read_csv('mcmaster.txt')

arrays=df.to_numpy()

objs=next(os.walk('.'))[1]

objs=objs[2:-1]

periods=[]

for j in objs:
    for i in arrays:
        if i[0]==j:
            periods.append(i[1])

    

#plot_fourier_grid(objs,periods,True)

obj=input('Select Cepheid: ')

period=0

for i in arrays:
    if i[0]==obj:
        period=i[1]

upper_bound=period*1.2
lower_bound=period/1.2

if obj in objs:
    print('~~~ '+obj+' ~~~\nLiterature Period Value: '+str(period))
    fourier_fitting(obj,period,2,5,True,False,20)
else:
    print('This object does not exist')



'''
RORY TO DO:

take median zpt error as min mag error

change folded axes to phase
color corrected vs raw mags
plot turning points as func of time
fix color plot bins and axes

biggest gaps in min and max fit vs data
then plot min/max difference on biggest diff in sampling in phase



MCMC

report notes:
period doubling - continual obs around peaks
1st rung on distance ladder
2-3 pages of intro/motivation
obs log appendix - separate legacy data
shift folding plots to minimum at origin
heliocentric corrections
photometry section in method
colour correction
light curve appendix
correlation coefficient for P-L relation
no error appendix - do within report

'''
