from imports import *
from astrolabs import get_data

def fourier_fitting(obj):
    times,mags,errors=get_data(obj)

    def fourier_function(t, *params):
        return params[0]+params[1]*np.sin(2*(np.pi/params[2])*t+params[3])+params[4]*np.sin(5*(np.pi/params[2])*t+params[5])
    
    fourier_values=[np.mean(mags),max(mags)-np.mean(mags),13.6,np.pi,max(mags)-np.mean(mags),np.pi]

    popt,cov=sp.optimize.curve_fit(fourier_function,times,mags,sigma=errors,p0=fourier_values,check_finite=True,maxfev=10**6)

    smooth_x=np.linspace(times[0], times[-1], 1000)

    output_popt=np.round(popt,2)

    plt.figure()
    plt.errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
    plt.plot(smooth_x,fourier_function(smooth_x, *popt),c='r',linestyle='dashed')
    plt.show()

    def chi_squared(model_params, model, x_data, y_data, y_err):
        return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
    reduced_chi=chi_squared(popt,fourier_function,times,mags,errors)/len(times)

    granularity=0.001
    error_period=popt[1]
    current_chi=reduced_chi

    while current_chi<reduced_chi+1:
        error_period+=granularity
        current_chi=chi_squared([popt[0],error_period,popt[1],popt[2],popt[3],popt[4],popt[5]],fourier_function,times,mags,errors)

    upper_error=abs(popt[1]-error_period)


    print('Period (days): '+str(output_popt[2])+' +/- '+str(upper_error))
    print('Fitted Function: '+str(output_popt[0])+'+'+str(output_popt[1])+'sin(2*2pi/'+str(output_popt[2])+'t+'+str(output_popt[3])+')+'+str(output_popt[4])+'sin(5*2pi/'+str(output_popt[2])+'t+'+str(output_popt[5])+')')
    print('Reduced Chi Squared: '+str(reduced_chi))

fourier_fitting('sz_cas')