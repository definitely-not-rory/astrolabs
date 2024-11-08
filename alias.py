from imports import *
from data_handling import get_data
from sinusoidal import sin_fitting
from fourier import fourier_fitting

def alias(obj,period,min_period,max_period,steps):

    times, mags, errors =get_data(obj)

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
    
    sin_bound_period=sin_fitting(obj,period)[0]
    fourier_bound_period=fourier_fitting(obj,period)[0]
    upper_bound=1.2*period
    lower_bound=period/1.2

    fig,ax=plt.subplots()
    plt.plot(periods,sin_chis,c='r')
    plt.plot(periods,fourier_chis,c='b')
    ax.axvline(sin_bound_period,c='r',linestyle='dashed')
    ax.axvline(sin_bound_period,c='b',linestyle='dashed')
    ax.axvline(upper_bound,c='k',linestyle='dashed')
    ax.axvline(lower_bound,c='k',linestyle='dashed')
    plt.xlabel('Period (days)')
    plt.ylabel('Chi Squared Value ( )')
    plt.show()

