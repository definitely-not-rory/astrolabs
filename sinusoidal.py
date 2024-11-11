from imports import *
from data_handling import get_data

def sin_fitting(obj,period):
    times, mags, errors =get_data(obj) #get data for given object

    amp=1
    period=period
    freq=2*np.pi/period
    phi=np.pi
    disp=np.mean(mags)

    initial_values = [amp,freq,phi,disp] #setting of trial values for both models

    def sin_function(t,*params):
        '''
        amplitude = params[0]
        frequency = params[1]
        phase = params[2]
        displacement = params[3]
        '''
        return params[0]*np.sin(params[1]*t-params[2])+params[3] #defining sine function for fitting
    
    amp_lo=-np.inf
    amp_hi=np.inf
    p_lo=period/1.2
    p_hi=period*1.2
    f_hi=2*np.pi/p_lo
    f_lo=2*np.pi/p_hi
    phi_lo=0
    phi_high=2*np.pi
    disp_lo=np.mean(mags)-1
    disp_hi=np.mean(mags)+1

    sin_bounds=([amp_lo,f_lo,phi_lo,disp_lo],[amp_hi,f_hi,phi_high,disp_hi])

    sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,bounds=sin_bounds,check_finite=True, maxfev=10**6)
    
    smooth_x=np.linspace(times[0], times[-1], 1000) #define x-range for plotting

    plt.errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3) #display raw data
    
    plt.plot(smooth_x,sin_function(smooth_x, *sinpopt),c='r',linestyle='dashed')

    plt.xlabel('Time (days)') #axes errors
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    plt.show()

    def chi_squared(model_params, model, x_data, y_data, y_err):
        return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
    sin_chi_val=chi_squared(sinpopt, sin_function, times, mags, errors)
    reduced_sin_chi=sin_chi_val/len(times)

    granularity=0.001
    error_freq=sinpopt[1]
    current_chi=reduced_sin_chi

    while current_chi<reduced_sin_chi+1:
        error_freq+=granularity
        current_chi=chi_squared([sinpopt[0],error_freq,sinpopt[1],sinpopt[2],sinpopt[3]],sin_function,times,mags,errors)

    upper_error=abs(sinpopt[1]-error_freq)

    error_freq=sinpopt[1]
    current_chi=reduced_sin_chi

    while current_chi<reduced_sin_chi+1:
        error_freq-=granularity
        current_chi=chi_squared([sinpopt[0],error_freq,sinpopt[1],sinpopt[2],sinpopt[3]],sin_function,times,mags,errors)

    lower_error=abs(sinpopt[1]-error_freq)

    mean_error=np.mean([lower_error,upper_error])

    fitted_period=2*np.pi/sinpopt[1]
    err_period=mean_error/sinpopt[1]**2

    print('\n~~~ '+obj+' Time Averaged Frequency and Period Data ~~~')
    print('\n~~~ Sinusoidal Model ~~~\nSin Frequency: '+str(sinpopt[1])+' +/- '+str(mean_error))
    print('Sin Period: '+str(fitted_period)+' +/- '+str(err_period))
    print('Sinusoidal Reduced Chi Squared: '+str(reduced_sin_chi)+'\n')

    folded_times=times
    for i in range(len(folded_times)):
        folded_times[i]=folded_times[i]%fitted_period
    
    print(folded_times)

    plt.figure()
    plt.errorbar(folded_times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
    plt.ylabel('Magnitude')
    plt.xlabel('Phase')
    plt.gca().invert_yaxis()
    plt.show()


    return fitted_period,err_period,reduced_sin_chi
