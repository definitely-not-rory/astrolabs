from imports import *

def snr_plot():
    apertures=np.array([0.8,1.8,2.8,3.8,4.8,5.8,6.8,7.8,10.8,12.8,15.8])
    signals=np.array([22232,140750,223410,262220,280470,289530,294630,297710,301710,303040,305210])
    skys=np.array([349.89,1823.6,4477.3,8194.7,13411,19924,27587,36458,70649,100350,153810])
    
    snrs=(signals-skys)/np.sqrt(signals+skys)

    plt.figure()
    plt.scatter(apertures,snrs,marker='x',c='k')
    plt.xlabel('Aperture Photometry Semi-Major Axis (px)')
    plt.ylabel('Signal-to-Noise Ratio ( )')
    plt.xlim((0,17))
    plt.ylim((0,520))
    plt.show()
