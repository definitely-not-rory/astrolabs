#IMPORTS
from imports import *

#DATA PROCESSING
def get_data_aphot(obj,use_mean):
    dates =os.listdir(obj+'/') #get list of available nights for given object
    times=[] #storage for data
    mags=[]
    errors=[]
    meantimes=[]
    meanmags=[]
    meanerrors=[]
    try: #try/except statement to ignore non-directory files/photometry files
        for date in dates: #iterate through available data
            df = pd.read_csv(obj+'/'+date+'/results.diff', delimiter=' ') #read results file outputted from raw2dif.py
            arrays=df.to_numpy()[:,:3] #remove NaN values (idk why they're there)
            meantimes=np.append(meantimes,np.mean(arrays[:,0]))
            meanmags=np.append(meanmags,np.median(arrays[:,1]))
            meanerrors=np.append(meanerrors,np.std(arrays[:,2])/len(arrays[:,2]))
            for point in arrays: #iterate through extracted data and store in correct places
                times.append(point[0])
                mags.append(point[1])
                errors.append(point[2])
    except:
        pass #ignore all non-directories
    times-=times[0] #normalise time values to begin at 0
    if use_mean==True:
        return meantimes,meanmags,meanerrors
    else:
        return times, mags, errors

#POST DATA COLLECTION CHANGES DATA PROCESSING
def get_data_manual(obj):
    df=pd.read_csv(obj+'/data.txt',delimiter=' ',header=None)
    arrays=df.to_numpy()
    times=arrays[:,0]
    counts=arrays[:,1]
    return times, counts

def calculate_magnitude(counts,zpt):
    mags=[]
    
    for i in range(len(counts)):
        mags=np.append(mags,zpt[i]-2.5*np.log10(counts[i]))
    return mags

def convert_to_times_mags(obj):
    times, counts = get_data_manual(obj)
    modified_times=[]
    for stringtime in times:
        objecttime=time.Time(stringtime, format='isot',scale='utc')
        modified_times.append(objecttime.mjd)
    zeros = []
    err_zeros = []
    with open(obj+'/zeros.txt','r') as readfile:
        for line in readfile:
            line = line.split()
            zeros = np.append(zeros,line[2])
            err_zeros = np.append(err_zeros,line[5])

    zeros=np.float64(zeros)
    err_zeros=np.float64(err_zeros)
    mags=calculate_magnitude(counts,zeros)
    errors=[]
    for i in range(len(counts)):
        errors=np.append(errors,np.sqrt(err_zeros[i]**2+(2.5*(np.sqrt(counts[i]))/(np.log(10)*counts[i]))**2))
    modified_times-=modified_times[0]
    return modified_times, mags, errors
    
#INITIAL PLOTS
def raw_plot(obj):
    fig, ax=plt.subplots()
    times, mags, errors =get_data_aphot(obj,False)
    markers,bars,caps=ax.errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=3)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    ax.set(xlabel='Time (days)',ylabel='Magnitude')
    fig.set_figheight(5)
    fig.set_figwidth(7.5)
    plt.title(obj)
    plt.show()

#CURVE FITTING
def fitting(obj,use_mean,newmethod):
    if newmethod==True:
        times,mags,errors=convert_to_times_mags(obj)
    else:
        times, mags, errors =get_data_aphot(obj,use_mean) #get data for given object

    initial_values = [max(mags)-(max(mags)+min(mags))/2,0.5,2*np.pi,(max(mags)+min(mags))/2] #setting of trial values for both models
    
    bounds = ([-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf])

    def sin_function(t,*params):
        '''
        amplitude = params[0]
        period = params[1]
        phase = params[2]
        displacement = params[3]
        '''
        return params[0]*np.sin(params[1]*t-params[2])+params[3] #defining sine function for fitting
   
    def sawtooth_function(t, *params):
        return params[0]*sp.signal.sawtooth(params[1]*t-params[2])+params[3] #defining sawtooth function for fitting
    

    sawpopt, sawcov = sp.optimize.curve_fit(sawtooth_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6, bounds=bounds) #run fitting for each model
    sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6, bounds=bounds)
    
    smooth_x=np.linspace(times[0], times[-1], 1000) #define x-range for plotting

    plt.errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3) #display raw data

    plt.plot(smooth_x,sawtooth_function(smooth_x, *sawpopt),c='r',linestyle='dashed') #plot fitted models
    plt.plot(smooth_x,sin_function(smooth_x, *sinpopt),c='b',linestyle='dashed')

    plt.xlabel('Time (days)') #axes errors
    plt.ylabel('Magnitude')

    sawpopt_errs = np.sqrt(np.diag(sawcov)) #calculate errors
    sinpopt_errs = np.sqrt(np.diag(sincov))

    def chi_squared(model_params, model, x_data, y_data, y_err):
        return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
    sin_chi_val=chi_squared(sinpopt, sin_function, times, mags, errors)
    reduced_sin_chi=sin_chi_val/len(times)

    saw_chi_val=chi_squared(sawpopt, sawtooth_function, times, mags, errors)
    reduced_saw_chi=saw_chi_val/len(times)

    raw_vs_average=''
    aphot_vs_hand=''

    if use_mean==True:
        raw_vs_average='Time Averaged'
    
    if newmethod==False:
        aphot_vs_hand='from aphot.py'

    print('\n~~~ '+raw_vs_average+'Frequency and Period Data'+aphot_vs_hand+' ~~~')
    print('\n ~~~ Sawtooth Model ~~~\nSawtooth Frequency: '+str(sawpopt[1]))
    print('Sawtooth Period: ' +str(2*np.pi/sawpopt[1])+ ' +/- '+str(sawpopt_errs[1]/sawpopt[1]**2)) #print calculated periods with errors
    print('Sawtooth Reduced Chi Squared: '+str(reduced_saw_chi)+'\n')

    print('\n~~~ Sinusoidal Model ~~~\nSin Frequency: '+str(sinpopt[1]))
    print('Sin Period: '+str(2*np.pi/sinpopt[1])+' +/- '+str(sinpopt_errs[1]/sinpopt[1]**2))
    print('Sinusoidal Reduced Chi Squared: '+str(reduced_sin_chi))

    plt.show()
