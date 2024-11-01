#IMPORTS
from imports import *

#DATA PROCESSING
def get_data(obj):
    dates =os.listdir(obj+'/') #get list of available nights for given object
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
   
    def sawtooth_function(t, *params):
        return params[0]*sp.signal.sawtooth(params[1]*t-params[2])+params[3] #defining sawtooth function for fitting
    
    def fourier(t, *params):
        pass

    def generate_initial_fourier(n):
        initial_fourier=np.ones(2*n+1)
        initial_fourier[0]=0.5
        return initial_fourier
    
    initial_fourier=generate_initial_fourier(1)

    def generate_fourier_bounds(input_params):
        lower_bound=-np.inf*np.ones(len(input_params))
        upper_bound=np.inf*np.ones(len(input_params))
        bounds=(lower_bound,upper_bound)
        return bounds

    sawpopt, sawcov = sp.optimize.curve_fit(sawtooth_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6) #run fitting for each model
    sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6)
    #fourierpopt, fouriercov = sp.optimize.curve_fit(fourier,times,mags,sigma=errors,absolute_sigma=True,p0=initial_fourier,check_finite=True,maxfev=10**6,bounds=generate_fourier_bounds(initial_fourier))
    
    smooth_x=np.linspace(times[0], times[-1], 1000) #define x-range for plotting

    plt.errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3) #display raw data
    
    plt.plot(smooth_x,sawtooth_function(smooth_x, *sawpopt),c='r',linestyle='dashed') #plot fitted models
    plt.plot(smooth_x,sin_function(smooth_x, *sinpopt),c='b',linestyle='dashed')
    #plt.plot(smooth_x,fourier(smooth_x,*fourierpopt),c='g',linestyle='dashed')

    plt.xlabel('Time (days)') #axes errors
    plt.ylabel('Magnitude')

    sawpopt_errs = np.sqrt(np.diag(sawcov)) #calculate errors
    sinpopt_errs = np.sqrt(np.diag(sincov))
    #fourierpopt_errs=np.sqrt(np.diag(fouriercov))

    def chi_squared(model_params, model, x_data, y_data, y_err):
        return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
    sin_chi_val=chi_squared(sinpopt, sin_function, times, mags, errors)
    reduced_sin_chi=sin_chi_val/len(times)

    saw_chi_val=chi_squared(sawpopt, sawtooth_function, times, mags, errors)
    reduced_saw_chi=saw_chi_val/len(times)

    #fourier_chi_val=chi_squared(fourierpopt,fourier,times,mags,errors)
    #reduced_fourier_chi=fourier_chi_val/len(times)

    raw_vs_average=''
    aphot_vs_hand=''

    print('\n~~~ Time Averaged Frequency and Period Data ~~~')
    print('\n ~~~ Sawtooth Model ~~~\nSawtooth Frequency: '+str(sawpopt[1]))
    print('Sawtooth Period: ' +str(2*np.pi/sawpopt[1])+ ' +/- '+str(sawpopt_errs[1]/sawpopt[1]**2)) #print calculated periods with errors
    print('Sawtooth Reduced Chi Squared: '+str(reduced_saw_chi)+'\n')

    print('\n~~~ Sinusoidal Model ~~~\nSin Frequency: '+str(sinpopt[1]))
    print('Sin Period: '+str(2*np.pi/sinpopt[1])+' +/- '+str(sinpopt_errs[1]/sinpopt[1]**2))
    print('Sinusoidal Reduced Chi Squared: '+str(reduced_sin_chi)+'\n')

    '''print('\n~~~ Fourier Model ~~~\nFourier Frequency: '+str(fourierpopt[0]))
    print('Fourier Period: '+str(2*np.pi/fourierpopt[0])+' +/- '+str(fourierpopt_errs[0]/fourierpopt[0]**2))
    print('Fourier Reduced Chi Squared: '+str(reduced_fourier_chi))'''

    plt.show()


'''

FIX CHI SQ

ALIASING PLOT - Chi vs Period, plot over number of observations

RESIDUALS

'''
