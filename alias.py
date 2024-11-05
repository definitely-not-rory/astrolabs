from imports import *
from data_handling import get_data

def alias(obj,min_period,max_period,steps):

    times, mags, errors =get_data(obj)

    periods=np.linspace(min_period,max_period,steps)
    sin_chis=[]
    fourier_chis=[]

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
        
        def fourier_function(t, *params):
            return params[0]+params[1]*np.sin(2*(np.pi/period)*t+params[2])+params[3]*np.sin(5*(np.pi/period)*t+params[4])
    
        fourier_values=[np.mean(mags),max(mags)-np.mean(mags),np.pi,max(mags)-np.mean(mags),np.pi]
                
        sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6)

        popt,cov=sp.optimize.curve_fit(fourier_function,times,mags,sigma=errors,p0=fourier_values,check_finite=True,maxfev=10**6)

        def chi_squared(model_params, model, x_data, y_data, y_err):
                return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
        sin_chi_val=chi_squared(sinpopt, sin_function, times, mags, errors)
        reduced_sin_chi=sin_chi_val/len(times)

        fourier_chi_val=chi_squared(popt,fourier_function,times,mags,errors)/len(times)

        sin_chis=np.append(sin_chis,reduced_sin_chi)
        fourier_chis=np.append(fourier_chis,fourier_chi_val)
    
    initial_values = [max(mags)-np.mean(mags),0.5,2*np.pi,np.mean(mags)] #setting of trial values for both models
    
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
    
    fourier_values=[np.mean(mags),max(mags)-np.mean(mags),13.6,np.pi,max(mags)-np.mean(mags),np.pi]

    sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6)
    fourierpopt, fouriercov = sp.optimize.curve_fit(fourier_function,times,mags,sigma=errors,absolute_sigma=True,p0=fourier_values,check_finite=True,maxfev=10**6)
    
    plt.figure()
    plt.plot(periods,sin_chis,c='r')
    plt.plot(periods,fourier_chis,c='b')
    plt.axvline(2*np.pi/0.5,c='k',linestyle='dashed')
    plt.text(2*np.pi/0.5+0.1,max(sin_chis)-3,'Initial Sinusoidal Fitting Starting Point',c='k')
    plt.axvline(2*np.pi/sinpopt[1],c='k',linestyle='dashed')
    plt.axvline(fourierpopt[2], c='k',linestyle='dashed')
    plt.text(2*np.pi/sinpopt[1]+0.2,max(sin_chis)-1,'Minimised Sinusoidal '+'$\chi^2$'+' period = '+str(np.round(2*np.pi/sinpopt[1],2)),c='k')
    plt.text(2*np.pi/fourierpopt[2]+0.2,max(sin_chis)-1,'Minimised Fourier '+'$\chi^2$'+' period = '+str(np.round(fourierpopt[2],2)),c='k')
    plt.xlabel('Period (days)')
    plt.ylabel('Chi Squared Value ( )')
    plt.show()

