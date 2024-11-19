from imports import *
from data_handling import get_data, raw_plot
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
#overfitted_modes(objs,periods)
#lit_observed_plot(objs,periods)
snr_plot()

obj=input('Select Cepheid: ')

period=0

for i in arrays:
    if i[0]==obj:
        period=i[1]
    
print('~~~ '+obj+' ~~~\nLiterature Period Value: '+str(period))

#chi_obs(obj,period,20)

'''
RORY TO DO:

fourier pages

change plots designs

residuals

report notes:
period doubling
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
