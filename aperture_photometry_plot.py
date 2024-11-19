from imports import *

def snr_plot():
    apertures=np.array([0.8,1.8,2.8,3.8,4.8,5.8,6.8,7.8,10.8,12.8,15.8])
    signals=np.array([22232,140750,223410,262220,280470,289530,294630,297710,301710,303040,305210])
    skys=np.array([88.443,37.003,108.40,653.11,925.24,544.48,1006.7,1195.8,501.97,1179.4,364.32])

    snrs=signals/np.sqrt(signals+skys)

    plt.figure()
    plt.scatter(apertures,snrs,marker='x',c='k')
    plt.xlabel('Aperture Photometry Semi-Major Axis (px)')
    plt.ylabel('Signal-to-Noise Ratio ( )')
    plt.show()
