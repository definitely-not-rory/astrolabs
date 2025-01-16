#IMPORTS
from imports import *
from min_error_calcs import *

#DATA PROCESSING
def get_data(obj):
    try:
        dates =os.listdir(obj+'/') #get list of available nights for given object
    except:
        sys.exit('This object is not in the directory')
    times=[]
    mjdtimes=[]
    mags=[]
    errors=[]
    meantimes=[] #storage for data
    meanmags=[]
    meanerrors=[]
    days=[]
    for date in dates: #iterate through available data
        folder=False
        try:
            df = pd.read_csv(obj+'/'+date+'/results.diff', delimiter=' ') #read results file outputted from raw2dif.py
            radec_df=pd.read_csv(obj+'/'+date+'/var_sky_position',delimiter=' ',header=None)
            folder=True
        except:
            try:
                df = pd.read_csv(obj+'/'+date+'/V/results.diff', delimiter=' ') #read results file outputted from raw2dif.py
                radec_df=pd.read_csv(obj+'/'+date+'/V/var_sky_position',delimiter=' ',header=None)
                folder=True
            except:
                pass
        if folder==True:
            arrays=df.to_numpy()[:,:3] #remove NaN values (idk why they're there)
            ra_dec=radec_df.to_numpy()[0]
            ra=ra_dec[0]
            dec=ra_dec[1]
            mjd_time=arrays[:,0]
            hjd_time=[]
            for mjd in mjd_time:
                hjd_time.append(pyasl.helio_jd(mjd,ra,dec))
            np.array(hjd_time)
            min_error=get_min_error()
            mjdtimes=np.append(mjdtimes,mjd_time)
            times=np.append(times,hjd_time)
            mags=np.append(mags,arrays[:,1])
            errors=np.append(errors,arrays[:,2])
            meantimes=np.append(meantimes,np.mean(arrays[:,0]))
            days=np.append(days,np.mean(arrays[:,0]))
            meanmags=np.append(meanmags,np.mean(arrays[:,1]))
            stdev=np.std(arrays[:,1])
            meaninerror=np.mean(arrays[:,2])
            if stdev>meaninerror:
                error=np.std(arrays[:,1])
            else:
                error=np.mean(arrays[:,2])
            if error<min_error:
                error=min_error
            meanerrors=np.append(meanerrors,error)
    meantimes-=meantimes[0]
    mjdtimes-=mjdtimes[0]
    times-=times[0]
    
    plt.figure()
    plt.errorbar(times,mags,errors,c='r',alpha=.5,capsize=3,ecolor='r',fmt='o',zorder=1,marker='x')
    plt.errorbar(meantimes,meanmags,meanerrors,c='b',capsize=3,fmt='o', marker=' ',ecolor='k',zorder=2)
    plt.scatter(meantimes,meanmags,c='b', marker='x',linewidth=2,zorder=3)
    plt.ylabel('Magnitude')
    plt.xlabel('Time (days)')
    plt.gca().invert_yaxis()
    plt.show()
    return meantimes,meanmags,meanerrors,days

def get_bv(obj):
    bands=os.listdir(obj+'/24_11_22/')
    try:
        for band in bands:
            df=pd.read_csv(obj+'/24_11_22/'+band+'/results.diff',delimiter=' ')
            arrays=df.to_numpy()[:,:3]
            if band=='B':
                bs=arrays[:,1]
            else:
                vs=arrays[:,1]
    except:
        print('Failed')
    return bs, vs
            


#INITIAL PLOTS
def raw_plot(obj):
    fig, ax=plt.subplots()
    times, mags, errors, days =get_data(obj)
    markers,bars,caps=ax.errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=3)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    ax.set(xlabel='Time (days)',ylabel='Magnitude')
    fig.set_figheight(5)
    fig.set_figwidth(7.5)
    plt.title(obj)
    plt.show()