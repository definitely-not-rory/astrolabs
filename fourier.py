from imports import *
from data_handling import get_data

def fourier_fitting(obj,period):
    times,mags,errors=get_data(obj)

    mean_mag=np.mean(mags)
    amp1=1
    period=period
    phi1=np.pi
    amp2=1
    phi2=np.pi

    def fourier_function(t, *params):
        return params[0]+params[1]*np.sin(2*(np.pi/params[2])*t+params[3])+params[4]*np.sin(5*(np.pi/params[2])*t+params[5])
    
    fourier_values=[mean_mag,amp1,period,phi1,amp2,phi2]

    amp1_lo=-np.inf
    amp1_hi=np.inf
    amp2_lo=-np.inf
    amp2_hi=np.inf
    p_lo=period/1.2
    p_hi=period*1.2
    phi1_lo=0
    phi1_hi=2*np.pi
    phi2_lo=0
    phi2_hi=2*np.pi
    disp_lo=np.mean(mags)-1
    disp_hi=np.mean(mags)+1

    fourier_bounds=([disp_lo,amp1_lo,p_lo,phi1_lo,amp2_lo,phi2_lo],[disp_hi,amp1_hi,p_hi,phi1_hi,amp2_hi,phi2_hi])

    popt,cov=sp.optimize.curve_fit(fourier_function,times,mags,sigma=errors,p0=fourier_values,bounds=fourier_bounds,check_finite=True,maxfev=10**6)

    smooth_x=np.linspace(times[0], times[-1], 1000)

    output_popt=np.round(popt,2)

    plt.figure()
    plt.errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
    plt.plot(smooth_x,fourier_function(smooth_x, *popt),c='r',linestyle='dashed')
    plt.xlabel('Time (days)') 
    plt.ylabel('Magnitude')
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

    error_period=popt[1]
    current_chi=reduced_chi

    while current_chi<reduced_chi+1:
        error_period-=granularity
        current_chi=chi_squared([popt[0],error_period,popt[1],popt[2],popt[3],popt[4],popt[5]],fourier_function,times,mags,errors)

    lower_error=abs(popt[1]-error_period)

    fitted_period=popt[1]
    mean_error=np.mean([lower_error,upper_error])

    print('Period (days): '+str(output_popt[2])+' +/- '+str(mean_error))
    print('Fitted Function: '+str(output_popt[0])+'+'+str(output_popt[1])+'sin(2*2pi/'+str(output_popt[2])+'t+'+str(output_popt[3])+')+'+str(output_popt[4])+'sin(5*2pi/'+str(output_popt[2])+'t+'+str(output_popt[5])+')')
    print('Reduced Chi Squared: '+str(reduced_chi))

    folded_times=times%fitted_period

    plt.figure()
    plt.errorbar(folded_times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
    plt.ylabel('Magnitude')
    plt.xlabel('Phase')
    plt.show()

    return fitted_period,mean_error,reduced_chi
