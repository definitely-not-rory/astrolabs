from imports import *

def snr_plot():
    apertures=np.array([0.8,1.8,2.8,3.8,4.8,5.8,6.8,7.8,10.8,12.8,15.8,25.1])
    signals=np.array([22232,140750,223410,262220,280470,289530,294630,297710,301710,303040,305210,338900])
    skys=np.array([349.89,1823.6,4477.3,8194.7,13411,19924,27587,36458,70649,100350,153810,194743])
    
    snrs=(signals-skys)/np.sqrt(signals+skys)

    print(f'Difference between 5 and 15: {snrs[4]/snrs[-2]}')
    print(f'Difference between 15 and 25: {snrs[-2]/snrs[-1]}')

    plt.figure()
    plt.scatter(apertures,snrs,marker='x',c='k')
    plt.xlabel('Aperture Photometry Semi-Major Axis (px)')
    plt.ylabel('Signal-to-Noise Ratio ( )')
    plt.fill_betweenx((0,520),0,5,facecolor='r',alpha=.3)
    plt.fill_betweenx((0,520),15,25,facecolor='b',alpha=.3)
    plt.axvline(5,c='r',linestyle='dashed')
    plt.axvline(15,c='b',linestyle='dashed')
    plt.axvline(25,c='b',linestyle='dashed')
    plt.xlim((0,26))
    plt.ylim((0,520))
    plt.show()

    '''was done using 4s exp or RW Casseiopeia'''
