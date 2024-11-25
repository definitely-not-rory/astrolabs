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
    days=[]
    for date in dates: #iterate through available data
        folder=False
        try:
            df = pd.read_csv(obj+'/'+date+'/results.diff', delimiter=' ') #read results file outputted from raw2dif.py
            folder=True
        except:
            try:
                df = pd.read_csv(obj+'/'+date+'/V/results.diff', delimiter=' ') #read results file outputted from raw2dif.py
                folder=True
            except:
                pass
        if folder==True:
            arrays=df.to_numpy()[:,:3] #remove NaN values (idk why they're there)
            meantimes=np.append(meantimes,np.mean(arrays[:,0]))
            days=np.append(days,np.mean(arrays[:,0]))
            meanmags=np.append(meanmags,np.mean(arrays[:,1]))
            mean_in_error=np.mean(arrays[:,2])
            std_in_data=np.std(arrays[:,1])
            if mean_in_error>std_in_data:
                meanerrors=np.append(meanerrors,np.mean(arrays[:,2]))
            else:    
                meanerrors=np.append(meanerrors,np.std(arrays[:,1]))
    meantimes-=meantimes[0]
    return meantimes,meanmags,meanerrors,days

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