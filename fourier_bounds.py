from imports import *
from fourier import fourier_fitting

def bounds_plot(obj,period):
    percentage_bounds=np.linspace(5,100,20)
    
    periods=[]
    errors=[]
    chis=[]


    for bound in percentage_bounds:
        fitted_period,error_from_jackknifing,reduced_chi,mean_mag, mean_mag_error=fourier_fitting(obj,period,2,5,False,False,bound)
        periods.append(fitted_period)
        errors.append(error_from_jackknifing)
        chis.append(reduced_chi)
    
    plt.errorbar(percentage_bounds,periods,yerr=errors,fmt='o')
    plt.show()


