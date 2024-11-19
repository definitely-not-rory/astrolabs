from imports import *
from fourier import fourier_fitting
from data_handling import get_data

def lit_observed_plot(objs,periods):
    observeds=[]
    observeds_errors=[]
    pos=0
    for obj in objs:
        period=periods[pos]
        observeds.append(fourier_fitting(obj,period,2,5,False,False,20)[0])
        observeds_errors.append(fourier_fitting(obj,period,2,5,False,False,20)[1])
        pos+=1
    vals=np.linspace(min(periods),max(periods),1000)
    plt.figure()
    plt.errorbar(periods,observeds,yerr=observeds_errors,c='k',marker='x',capsize=3,fmt='o')
    plt.plot(vals,vals,c='r',linestyle='dashed')
    plt.show()
