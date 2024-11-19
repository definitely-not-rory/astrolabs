from imports import *
from fourier import fourier_fitting

def bounds_plot(obj,period):
    percentage_bounds=np.linspace(5,500,100)
    
    periods=[]
    errors=[]
    chis=[]


    for bound in percentage_bounds:
        fitted_period,error_from_jackknifing,reduced_chi,mean_mag, mean_mag_error=fourier_fitting(obj,period,2,5,False,False,bound)
        periods.append(fitted_period)
        errors.append(error_from_jackknifing)
        chis.append(reduced_chi)
    
    plt.errorbar(percentage_bounds,periods,yerr=errors,fmt='o',c='k',marker='x',capsize=3)
    plt.xlabel('Percentage of literature period bound used in $\chi^2$ minimisation fitting')
    plt.ylabel('Fitted Period from $\chi^2$ minimisation fitting')
    plt.axhline(period,c='r',linestyle='dashed')
    plt.show()


