#POST DATA COLLECTION CHANGES DATA PROCESSING
def get_data_manual(obj):
    df=pd.read_csv(obj+'/data.txt',delimiter=' ',header=None)
    arrays=df.to_numpy()
    times=arrays[:,0]
    counts=arrays[:,1]
    try:
        exposures=arrays[:,2]
        return times, counts, exposures
    except:
        return times, counts
    
def calculate_magnitude(counts,zpt):
    mags=[]
    for i in range(len(counts)):
        mags=np.append(mags,zpt[i]-2.5*np.log10(counts[i]))
    return mags

def convert_to_times_mags(obj):
    try:
        times, counts, exposures = get_data_manual(obj)
    except:
        times, counts = get_data_manual(obj)
    modified_times=[]
    for stringtime in times:
        objecttime=astrotime.Time(stringtime, format='isot',scale='utc')
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

