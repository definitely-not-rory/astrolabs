def get_automag_pars():
    file = '/mnt/64bin/auto_astrom/automag_driver'

    if not os.path.exists('automag_driver'):
        printlog('No automag_driver in current directory.  Using default at /mnt/64bin/auto_astrom/automag_driver')
        if not os.path.exists(file):
            printlog('No default autmag_driver file.  Will not continue.  Aborting...')
            sys.exit()

    if os.path.exists('automag_driver'):
        file = 'automag_driver'

    printlog('automag_driver: using file '+file)
    data = open(file,'r')
    n = len(open(file).readlines(  ))
    for i in range(0,n):
        line = data.readline()
        word = line.split()
        if word[0] == 'radius':
            rad = float(word[2])
        if word[0] == 'radius_1':
            r_in = float(word[2])
        if word[0] == 'radius_2':
            r_out = float(word[2])
    printlog('aperture parameters; rad, r_in, r_out: '+str(rad)+', '+str(r_in)+', '+str(r_out))
    return rad,r_in,r_out

# --------------------------------------------------------------------------------------------------------------

def get_minor_position(datetime):
    # print("getMJD(datetime) = ",getMJD(datetime))
    #timeStart=time.time()
    # Replace this subourine with the content of the actual file
    #ra_deg, dec_deg = subprocess.getoutput("python /mnt/64bin/get_jpl_position.py3 " + getMJD(datetime)).split()
    ra_deg, dec_deg = get_jpl_position(getMJD(datetime))
    
    ra_deg=float(ra_deg)
    dec_deg=float(dec_deg)
    printlog('RA/Dec of minor body: '+str(ra_deg)+' '+str(dec_deg))
    return ra_deg, dec_deg

# --------------------------------------------------------------------------------------------------------------

def get_jpl_position(MJD):
    req_jd  = float(MJD)
    """This is just '/mnt/64bin/get_jpl_position.py3'
        This is the half that is called every subroutine"""
    if req_jd < min(jtime) or req_jd > max(jtime):
        print("Requested MJD outside range in file: jpl_ephemeris")
        print("ABORTING")
        sys.exit()

    ihit = 0
    dmin = 99999999999.0
    for i in range(len(jtime)):
        if abs(jtime[i] - req_jd) < dmin:
            dmin = abs(jtime[i] - req_jd)
            ihit = i

    if ihit < 2:
        ihit = 2
    if ihit > len(jtime) - 3:
        ihit = len(jtime) - 3
    jjt = np.array(
        [jtime[ihit - 2], jtime[ihit - 1], jtime[ihit], jtime[ihit + 1], jtime[ihit + 2]]
    )
    jjra = np.array(
        [jra[ihit - 2], jra[ihit - 1], jra[ihit], jra[ihit + 1], jra[ihit + 2]]
    )
    jjdec = np.array(
        [jdec[ihit - 2], jdec[ihit - 1], jdec[ihit], jdec[ihit + 1], jdec[ihit + 2]]
    )
    jpl_ra = np.interp(req_jd, jjt, jjra)
    jpl_dec = np.interp(req_jd, jjt, jjdec)
    
    #recordRaDec(MJD,jpl_ra, jpl_dec) #!!! Added to compare Ra&Decs
    return jpl_ra, jpl_dec

# --------------------------------------------------------------------------------------------------------------

def checkargs():
    hostname = socket.gethostname()

    if len(sys.argv) != 5:
        print(" ")
        print("Input information required on the command line is")
        print("                  TELESCOPE  OBS_DAT  START  END")
        print(" ")
        print("  where TELESCOPE is 'draco2', 'east-14', 'west-14', etc")
        print("  OBS_DATE is the observation date given as 'yy_mm_dd'")
        print("  START is the first image number to process")
        print("  END is the last image number to process")
        print("  e.g.")
        print("    python3 /mnt/64bin/auto_astrom/aphot.py draco2 15_10_24 400 412")
        print(" ")
        sys.exit()

    TELESCOPE = sys.argv[1]
    if not TELESCOPE in [
            "draco2",
            "east-14",
            "east-16",
            "far-east-14",
            "far-east-16",
            "far-east-20",
            "west-12",
            "west-14",
            "14",
            "draco",
            "pt5m",
            "10_2",
            "12"]:
        print("Possible telescope names are:")
        print("          draco2  east-14  far-east-14  far-east-16 far-east-20 west-12  west-14 draco 14 pt5m")
        sys.exit()

    ODATE = sys.argv[2]
    if len(ODATE) != 8:
        print(" ")
        print("Date format required is yy_mm_dd, e.g.  14_11_05")
        print(" ")
        sys.exit()

    if int(ODATE[0:2]) > 80:
        OYEAR = "19" + ODATE[0:2]
    else:
        OYEAR = "20" + ODATE[0:2]
    if hostname == "capella.astrolab" or hostname == "canopus.astrolab":
        DATA_PATH = "/archive/" + TELESCOPE + "/" + OYEAR + "/" + ODATE + "/"
    else:
        DATA_PATH = "/mnt//archive/" + TELESCOPE + "/" + OYEAR + "/" + ODATE + "/"

    if not os.path.exists(DATA_PATH):
        print(" ")
        print("ERROR:  ", DATA_PATH, "does not exists")
        print(" ")
        sys.exit()

    return DATA_PATH

# --------------------------------------------------------------------------------------------------------------

def printlog(word):
    print(word)
    log.write(str(word)+'\n')
    return


# --------------------------------------------------------------------------------------------------------------

def readcoords_var(datetime,moving_target):

    if moving_target==1:
        #ra,dec = get_minor_position(datetime)
        datetime_MJD = getMJD(datetime)
        ra,dec = get_jpl_position(datetime_MJD)
        return ra,dec
    else:
        vfile = "var_sky_position"
        f = open(vfile,'r')
        a = f.readline()
        a = a.replace("\n",'')
        a = a.split(" ")
        f.close()
        ra = float(a[0])
        dec = float(a[1])
        return ra,dec
    
# --------------------------------------------------------------------------------------------------------------

def readcoords_cal(file):
    f = open(file,'r')
    a = f.readline()
    a = a.replace("\n",'')
    a = a.split(" ")
    b = f.readline()
    b = b.replace("\n",'')
    b = b.split(" ")
    f.close()
    return ((float(a[0]),float(a[1])),(float(b[0]),float(b[1])))

# --------------------------------------------------------------------------------------------------------------

def checkforfiles():
    cur_dir = os.getcwd()
    nspace = cur_dir.find(" ")
    if nspace > 0:
        print(" ")
        print("       The working directory has a space character.")
        print(" ")
        sys.exit()
        
    if not os.path.exists("cal_sky_positions"):
        print(" ")
        print("       Missing file: cal_sky_positions")
        print(" ")
        sys.exit()

    if not (os.path.exists("var_sky_position") or os.path.exists("jpl_ephemeris")):
        print(" ")
        print("       Missing file:")
        print("              Need either the 'var_sky_position' or 'jpl_ephemeris' file")
        print(" ")
        print("       For moving targets, create a 'jpl_ephemeris' file with:")
        print("       python /mnt/64bin/get_jpl_eph.py")
        sys.exit()

# --------------------------------------------------------------------------------------------------------------
        
def getscale(hdr):
    w = wcs.WCS(hdr)
    ra0,dec0 = w.wcs_pix2world(0,0,1)
    ra1,dec1 = w.wcs_pix2world(0,1,1)
    dra = (ra1-ra0)*np.cos(dec0*np.pi/180)
    ddec = (dec1-dec0)
    asecpix = (dra**2+ddec**2)**0.5
    #asecpix = ((dec1-dec0)**2+(ra0-ra1)**2)**0.5
    asecpix = abs(3600*asecpix)
    return asecpix

# ------------------------------------------------------------------------------------------------------------------------------------------

def fits_read(file,i):
    hdu = fits.open(file)
    data = hdu[i].data
    hdr = hdu[i].header
    return data,hdr

# ------------------------------------------------------------------------------------------------------------------------------------------

def ad2xy(hdr,ra,dec):
    w = wcs.WCS(hdr)
    xc,yc = w.wcs_world2pix(ra,dec,1)
    return xc,yc

# ------------------------------------------------------------------------------------------------------------------------------------------

def getMJD(datetime):
    t = Time(datetime, format='isot', scale='utc')
    date_mjd = "{:.5f}".format(t.mjd)
    return(date_mjd)

# ------------------------------------------------------------------------------------------------------------------------------------------

def checkxy(x,y,sz):
    flag = 0
    if x <= 0 or x > sz[1] or y <= 0 or y >= sz[0] or np.isfinite(x)==False or np.isfinite(y)==False:
        flag = 1
        x = 0
        y = 0
        printlog('WARNING: RA/Dec are not in image range')
    return x,y,flag

# -----------------------------------------------------------------------------------------------------------------------------------------

def getphot(img,hdr,radecvar,radec_cal1,radec_cal2,rad_asec,r_in_asec,r_out_asec,file):

    # get pixel scale of image
    spix = getscale(hdr)
    printlog('pixel scale is: '+str(spix))

    # convert RA/Dec of variable star and two calibration stars to x/y coordinates on this image
    x_var,y_var   = ad2xy(hdr,radec_var[0],radec_var[1])
    x_cal1,y_cal1 = ad2xy(hdr,radec_cal1[0],radec_cal1[1])
    x_cal2,y_cal2 = ad2xy(hdr,radec_cal2[0],radec_cal2[1])
    printlog('x/y of var, cal1, cal2: '+str(x_var)+' '+str(y_var)+' '+str(x_cal1)+' '+str(y_cal1)+' '+str(x_cal2)+' '+str(y_cal2))

    # convert aperture sizes from arcseconds to pixels
    rad = rad_asec / spix        # pixels
    r_in = r_in_asec / spix      # pixels
    r_out = r_out_asec / spix    # pixels
    
    # CCD parameters (used for noise calculation)
    gain = 1.
    dark = 0
    rd = 15.

    # check if x/y coordinates are on the frame.
    sz = img.shape
    x_var,y_var,xy_var_bad = checkxy(x_var,y_var,sz)
    x_cal1,y_cal1,xy_cal1_bad = checkxy(x_cal1,y_cal1,sz)
    x_cal2,y_cal2,xy_cal2_bad = checkxy(x_cal2,y_cal2,sz)

    # if any of the x/y positions of the variable or two calib stars are not in the image, no photometry is calculated.
    abort = 0
    if xy_var_bad == 1 or xy_cal1_bad == 1 or xy_cal2_bad == 1:
        word = file+' ' + str(x_var) +' '+ str(y_var)+' '+str(x_cal1)+' '+str(y_cal1)+' '+str(x_cal2)+' '+str(y_cal2)
        g.write(word+'\n')
        abort=1
        return 0,abort

    if abort == 0:
        # create tuple of x/ypositions
        positions = [(x_var,y_var),(x_cal1,y_cal1),(x_cal2,y_cal2)]
        # create circular aperture (for photometry calculation)
        aperture = CircularAperture(positions, r=rad)
        # create annulus aperture (for sky calculation)
        annulus_aperture = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
        # create mask
        annulus_masks = annulus_aperture.to_mask(method='center')

        # calculate the data and (sigma clipped) sky
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(img)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        # write results to output object
        phot = aperture_photometry(img, aperture)
        phot['annulus_median'] = bkg_median
        phot['aper_bkg'] = bkg_median * aperture.area
        phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
        # noise calculation
        phot['aper_sum_bkgsub_err'] = np.sqrt((phot['aperture_sum']) + (phot['aper_bkg']) + dark + rd**2)
        for col in phot.colnames:
            phot[col].info.format = '%.8g'  # for consistent table output
        printlog(phot)
        return(phot,abort)

# ------------------------------------------------------------------------------------------------------------------------------------------

import socket, math
from ephem import EllipticalBody, Observer
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import wcs
import photutils
from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.datasets import make_100gaussians_image
#from photutils import aperture_photometry
#from photutils import CircularAperture, CircularAnnulus
#from photutils.datasets import make_100gaussians_image
import numpy as np
import sys
from astropy.time import Time
import os.path
import os, subprocess
from scipy.interpolate import interp1d

checkforfiles()
data_path = checkargs()

log = open('aphot.log','w')

#######################

flatfile = 'master_flat.fits'
# check for flatfield
flat_exist = os.path.isfile(flatfile)
if flat_exist == 1:
    flat,flathdr = fits_read(flatfile,0)

# read calibration star RA/Dec positions
radec_cal1,radec_cal2 = readcoords_cal('cal_sky_positions')
printlog('read ra/dec for calibration star of: '+str(radec_cal1)+' '+str(radec_cal2))

rad,r_in,r_out = get_automag_pars()
#rad = 5.    # arcsec, aperture
#r_in = 15   # arcsec, sky inner annulus
#r_out = 25. # arcsec, sky outer annulus
#sys.exit()
# random zero-point constant
zpt = 20.0

#######################

moving_target = 0
if os.path.exists("jpl_ephemeris"):
    printlog("found jpl_ephemeris")
    moving_target=1
    # read jpl_ephemeris file
    jtime = []
    jra = []
    jdec = []
    f = open("jpl_ephemeris", "r")
    line = f.readline()
    while True:
        line = f.readline()
        if line == "":
            break
        word = line.split()
        jtime.append(float(word[0]))
        jra.append(float(word[2]))
        jdec.append(float(word[3]))
    f.close()
    j_start = min(jtime)
    j_end = max(jtime)



f = open('summary.obs','a')
g = open('FAILURES.obs','w')

istart = int(sys.argv[3])
iend = int(sys.argv[4])

for i in range(istart,iend):
    #if i <= 9: pad = '000'
    #if i >= 10 and i <= 99: pad = '00'
    #if i >= 100 and i <= 999: pad = '0'
    #if i >= 1000: pad = ''
    printlog("-----------------------------------------")
    printlog(i)

    #file = 'data/ad' + pad + str(i) + '.fits'
    # data_path='data'
    file = data_path+'/adbr0' + str(i) + '.fits'
    print(file)
    #sys.exit()
    
    # check for file and process if it exists
    exist = os.path.isfile(file)
    if exist == 1:
        printlog(file)
        img,hdr = fits_read(file,0)
        datetime = hdr['DATE-OBS']
    
        # divide by flat
        if flat_exist == 1:
            img = img/flat
            printlog('flat found.  dividing by flatfield')
        if flat_exist == 0:
            printlog('no flat found.  will not flatfield')

        # read variable star position
        radec_var = readcoords_var(datetime,moving_target)
        printlog('read ra/dec for target of: '+str(radec_var))

        # calculate photometry and also get "flag".  if flag=0, photometry was successful.  if flag=1, an error occured.
        phot,flag = getphot(img,hdr,radec_var,radec_cal1,radec_cal2,rad,r_in,r_out,file)

        # if photometry was successful, write results to summary.obs
        if flag==0:
            MJD = getMJD(datetime)
            #printlog(MJD)

            mag_var = zpt-2.5*np.log10(phot[0]['aper_sum_bkgsub'])
            mag_var_err = zpt-2.5*np.log10(phot[0]['aper_sum_bkgsub']-phot[0]['aper_sum_bkgsub_err'])-mag_var

            mag_cal1 = zpt-2.5*np.log10(phot[1]['aper_sum_bkgsub'])
            mag_cal1_err = zpt-2.5*np.log10(phot[1]['aper_sum_bkgsub']-phot[1]['aper_sum_bkgsub_err'])-mag_cal1

            mag_cal2 = zpt-2.5*np.log10(phot[2]['aper_sum_bkgsub'])
            mag_cal2_err = zpt-2.5*np.log10(phot[2]['aper_sum_bkgsub']-phot[2]['aper_sum_bkgsub_err'])-mag_cal2
    
            word = str(MJD)+' '+"{:.4f}".format(mag_var)+' '+"{:.4f}".format(mag_var_err)+\
                   ' '+"{:.4f}".format(mag_cal1)+' '+"{:.4f}".format(mag_cal1_err)+' '+\
                   "{:.4f}".format(mag_cal2)+' '+"{:.4f}".format(mag_cal2_err)+' * ' + \
                   file + ' ' + datetime[0:10]+' '+datetime[11:25]
            f.write(word+'\n')
f.close
g.close
log.close

# Changes made October 2021:
# np.float replaced with float
# Modified: get_minor_position & readcoords_var
