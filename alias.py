from imports import *
from data_handling import get_data
from sinusoidal import sin_fitting
from fourier import fourier_fitting

def alias(obj,period,min_period,max_period,steps,show_plots):

    times, mags, errors,days =get_data(obj)

    def sin_function(t,*params):
            '''
            amplitude = params[0]
            frequency = params[1]
            phase = params[2]
            displacement = params[3]
            '''
            return params[0]*np.sin(params[1]*t-params[2])+params[3] #defining sine function for fitting
    
    def fourier_function(t, *params):
        return params[0]+params[1]*np.sin(2*(np.pi/params[2])*t+params[3])+params[4]*np.sin(5*(np.pi/params[2])*t+params[5])
    
    def chi_squared(model_params, model, x_data, y_data, y_err):
                return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)

    periods=np.linspace(min_period,max_period,steps)
    sin_chis=[]
    fourier_chis=[]

    for a_period in periods:
        a_amp=1
        a_period=a_period
        a_freq=2*np.pi/a_period
        a_phi=np.pi
        a_disp=np.mean(mags)

        def alias_sin_function(t,*params):
            '''
            amplitude = params[0]
            frequency = params[1]
            phase = params[2]
            displacement = params[3]
            '''
            return params[0]*np.sin(a_freq*t-params[1])+params[2] #defining sine function for fitting

        a_initial_values = [a_amp,a_phi,a_disp] #setting of trial values for both models

        a_sinpopt, a_sincov = sp.optimize.curve_fit(alias_sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=a_initial_values,check_finite=True, maxfev=10**6)

        a_mean_mag=np.mean(mags)
        a_amp1=1
        a_period=a_period
        a_phi1=np.pi
        a_amp2=1
        a_phi2=np.pi

        def alias_fourier_function(t, *params):
            return params[0]+params[1]*np.sin(2*(np.pi/a_period)*t+params[2])+params[3]*np.sin(5*(np.pi/period)*t+params[4])

        a_fourier_values=[a_mean_mag,a_amp1,a_phi1,a_amp2,a_phi2]

        a_fourier_popt,a_fourier_cov=sp.optimize.curve_fit(alias_fourier_function,times,mags,sigma=errors,p0=a_fourier_values,check_finite=True,maxfev=10**6)
    
        a_sin_chi_val=chi_squared(a_sinpopt, alias_sin_function, times, mags, errors)
        a_reduced_sin_chi=a_sin_chi_val/len(times)

        a_fourier_chi_val=chi_squared(a_fourier_popt,alias_fourier_function,times,mags,errors)/len(times)

        sin_chis=np.append(sin_chis,a_reduced_sin_chi)
        fourier_chis=np.append(fourier_chis,a_fourier_chi_val)
    
    fourier_bound_period=fourier_fitting(obj,period,2,5,False,False,20)[0]
    bound_error=fourier_fitting(obj,period,2,5,False,False,20)[1]
    upper_bound=1.2*period
    lower_bound=period/1.2
    delta_chis=fourier_chis-min(fourier_chis)

    min_index=np.argmin(delta_chis)
    '''
    low_pos=min_index
    found_low_err=False
    while found_low_err==False:
        if delta_chis[low_pos]>1:
              found_low_err=True
        else:
              low_pos-=1
    
    high_pos=min_index
    found_high_err=False
    while found_high_err==False:
        if delta_chis[high_pos]>1:
              found_high_err=True
        else:
              high_pos+=1
    
    low_err=periods[min_index]-periods[low_pos]
    high_err=periods[high_pos]-periods[min_index]

    print(low_err)
    print(high_err)'''

    if show_plots==True:
        fig,ax=plt.subplots()
        plt.plot(periods,delta_chis,c='r')
        ax.axvline(fourier_bound_period,c='b',linestyle='dashed')
        ax.axvline(fourier_bound_period-bound_error,c='b',linestyle='dashed',alpha=0.5)
        ax.axvline(fourier_bound_period+bound_error,c='b',linestyle='dashed',alpha=0.5)
        ax.axhline(1,c='g',linestyle='dashed')
        ax.axvline(min_period,c='k',linestyle='dashed')
        ax.axvline(max_period,c='k',linestyle='dashed')
        #ax.axvline(fourier_bound_period-low_err,c='g',linestyle='dashed',alpha=.5)
        #ax.axvline(fourier_bound_period+high_err,c='g',linestyle='dashed',alpha=.5)
        plt.ylabel('$\Delta\chi^2$')
        plt.xlabel('Period (days)')
        plt.show()
    #return low_err, high_err


