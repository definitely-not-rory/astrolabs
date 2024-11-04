#IMPORTS
from imports import *

#DATA PROCESSING
def get_data(obj):
    try:
        dates =os.listdir(obj+'/') #get list of available nights for given object
    except:
        sys.exit('This object is not in the directory')
    meantimes=[] #storage for data
    meanmags=[]
    meanerrors=[]
    for date in dates: #iterate through available data
        folder=False
        try:
            df = pd.read_csv(obj+'/'+date+'/results.diff', delimiter=' ') #read results file outputted from raw2dif.py
            folder=True
        except:
            pass
        if folder==True:
            arrays=df.to_numpy()[:,:3] #remove NaN values (idk why they're there)
            meantimes=np.append(meantimes,np.mean(arrays[:,0]))
            meanmags=np.append(meanmags,np.mean(arrays[:,1]))
            mean_in_error=np.mean(arrays[:,2])
            std_in_data=np.std(arrays[:,1])
            if mean_in_error>std_in_data:
                meanerrors=np.append(meanerrors,np.mean(arrays[:,2]))
            else:    
                meanerrors=np.append(meanerrors,np.std(arrays[:,1]))
    meantimes-=meantimes[0]
    return meantimes,meanmags,meanerrors

#INITIAL PLOTS
def raw_plot(obj):
    fig, ax=plt.subplots()
    times, mags, errors =get_data(obj)
    markers,bars,caps=ax.errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=3)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    ax.set(xlabel='Time (days)',ylabel='Magnitude')
    fig.set_figheight(5)
    fig.set_figwidth(7.5)
    plt.title(obj)
    plt.show()

#CURVE FITTING
def fitting(obj):
    times, mags, errors =get_data(obj) #get data for given object
    initial_values = [max(mags)-np.mean(mags),0.5,2*np.pi,np.mean(mags)] #setting of trial values for both models

    def sin_function(t,*params):
        '''
        amplitude = params[0]
        frequency = params[1]
        phase = params[2]
        displacement = params[3]
        '''
        return params[0]*np.sin(params[1]*t-params[2])+params[3] #defining sine function for fitting
   
    def generate_initial_fourier(n):
        initial_fourier=np.ones(2*n+1)
        initial_fourier[0]=0.5
        return initial_fourier
    
    sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6)
    
    smooth_x=np.linspace(times[0], times[-1], 1000) #define x-range for plotting

    plt.errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3) #display raw data
    
    plt.plot(smooth_x,sin_function(smooth_x, *sinpopt),c='b',linestyle='dashed')

    plt.xlabel('Time (days)') #axes errors
    plt.ylabel('Magnitude')
    plt.show()
    sinpopt_errs = np.sqrt(np.diag(sincov))

    def chi_squared(model_params, model, x_data, y_data, y_err):
        return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
    sin_chi_val=chi_squared(sinpopt, sin_function, times, mags, errors)
    reduced_sin_chi=sin_chi_val/len(times)

    granularities=np.linspace(0.001,0.01,10)

    for freq_granularity in granularities:

        error_reduced_sin_chi_val=reduced_sin_chi
        error_freq=sinpopt[1]
        upper_freqs=[error_freq]
        reduced_chis_upper=[error_reduced_sin_chi_val]

        while error_reduced_sin_chi_val<reduced_sin_chi+1:
            error_freq+=freq_granularity
            upper_freqs.append(error_freq)
            error_reduced_sin_chi_val=chi_squared([sinpopt[0],error_freq,sinpopt[2],sinpopt[3]],sin_function,times,mags,errors)/len(times)
            reduced_chis_upper.append(error_reduced_sin_chi_val)
        
        upper_error=error_freq-sinpopt[1]

        error_reduced_sin_chi_val=reduced_sin_chi
        error_freq=sinpopt[1]
        lower_freqs=[error_freq]
        reduced_chis_lower=[error_reduced_sin_chi_val]

        while error_reduced_sin_chi_val<reduced_sin_chi+1:
            error_freq-=freq_granularity
            lower_freqs.append(error_freq)
            error_reduced_sin_chi_val=chi_squared([sinpopt[0],error_freq,sinpopt[2],sinpopt[3]],sin_function,times,mags,errors)/len(times)
            reduced_chis_lower.append(error_reduced_sin_chi_val)

        lower_error=sinpopt[1]-error_freq

        mean_error_in_freq=np.mean([upper_error,lower_error])    

        plt.plot(upper_freqs,reduced_chis_upper,c='r',alpha=1.1-(freq_granularity*100))
        plt.plot(lower_freqs,reduced_chis_lower,c='b',alpha=1.1-(freq_granularity*100))
    
    plt.axvline(sinpopt[1],c='k',linestyle='dashed')
    plt.axhline(reduced_sin_chi,c='k',linestyle='dashed')
    plt.axhline(reduced_sin_chi+1,c='k',linestyle='dashed')
    plt.xlabel('Frequency (1/days)')
    plt.ylabel('Chi Squared ( )')
    plt.show()

    #plt.figure()

    print('\n~~~ Time Averaged Frequency and Period Data ~~~')
    print('\n~~~ Sinusoidal Model ~~~\nSin Frequency: '+str(sinpopt[1])+' +/- '+str(mean_error_in_freq))
    print('Sin Period: '+str(2*np.pi/sinpopt[1])+' +/- '+str(upper_error/sinpopt[1]**2))
    print('Sinusoidal Reduced Chi Squared: '+str(reduced_sin_chi)+'\n')


fitting('sw_cas')

'''

ALIASING PLOT - Chi vs Period, plot over number of observations

RESIDUALS

'''
