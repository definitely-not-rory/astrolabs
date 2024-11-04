from imports import *
from astrolabs import get_data, fitting

def alias(obj,min_period,max_period,steps):

    fitting(obj)

    times, mags, errors =get_data(obj)

    periods=np.linspace(min_period,max_period,steps)
    chis=[]

    for period in periods:
        initial_values = [max(mags)-np.mean(mags),2*np.pi,np.mean(mags)] #setting of trial values for both models

        def sin_function(t,*params):
            '''
            amplitude = params[0]
            frequency = params[1]
            phase = params[2]
            displacement = params[3]
            '''
            return params[0]*np.sin((2*np.pi/period)*t-params[1])+params[2] #defining sine function for fitting
                
        sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6)

        def chi_squared(model_params, model, x_data, y_data, y_err):
                return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
        sin_chi_val=chi_squared(sinpopt, sin_function, times, mags, errors)
        reduced_sin_chi=sin_chi_val/len(times)

        chis=np.append(chis,reduced_sin_chi)

    plt.figure()
    plt.plot(periods,chis)
    plt.xlabel('Period (days)')
    plt.ylabel('Chi Squared Value ( )')
    plt.show()

alias('sz_cas',1,15,3000)

        