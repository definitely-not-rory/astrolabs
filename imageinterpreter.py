import matplotlib
from matplotlib import pyplot as plt
import astropy
from astropy.io import fits
import numpy as np

plt.close('all')
ff = 0.9

file = '24_10_17_sz_cas.fits'

hdu = fits.open(file)
img = hdu[0].data
hdr = hdu[0].header

mm = np.median(img)
gd = np.where(img)
allval = img[gd]
sort_allval = np.sort(allval)
n = len(sort_allval)
vhi = sort_allval[int(n*ff)]

plt.axis('off')

ii = plt.imshow(img,vmin=mm,vmax=vhi,origin='lower',cmap='Greys')
plt.savefig('24_10_17_sz_cas.png', bbox_inches='tight')
plt.show(block=True)

