from astrolabs import get_data, raw_plot
from imports import *

def fittingtesting():
    times, mags, errors = get_data('cg_cas',False)
    times, mags, errors = times[2:-7], mags[2:-7], errors[2:-7]

    plt.figure()
    plt.scatter(times, mags)
    initial_values=[1,1]
    def linear(x, *params):
        return params[0]*x+params[1]
    popt, cov = sp.optimize.curve_fit(linear,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True)

    smooth_x=np.linspace(times[0], times[-1], 1000) #define x-range for plotting

    plt.plot(smooth_x,linear(smooth_x,*popt))
    print('gradient: '+str(popt[0]))
    print('intercept: '+str(popt[1]))

    fitted_params=[popt[0],popt[1]]
    print(fitted_params)

def aphot_cg_cas_testing():
    times, mags, errors = get_data('cg_cas',False)
    hand_values= [11.758,11.725]
    hand_errors=[0.019,0.018]
    differences=[abs(hand_values[0]-mags[0]),abs(hand_values[1]-mags[1])]
    error_differences=[abs(hand_errors[0]-errors[0]),abs(hand_errors[1]-errors[1])]
    difference_in_points=[abs(mags[0]-mags[1]),abs(hand_values[0]-hand_values[1])]
    mag_ratios=[mags[0]/hand_values[0],mags[1]/hand_values[1]]
    err_ratios=[errors[0]/hand_errors[0],errors[1]/hand_errors[1]]

    print('\n~~~~~ ad0101.fits ~~~~~\naphot: '+str(mags[0])+' +/- '+str(errors[0])+'\nby hand: '+str(hand_values[0])+' +/- '+str(hand_errors[0]))
    print('\n~~~~~ ad0102.fits ~~~~~\naphot: '+str(mags[1])+' +/- '+str(errors[1])+'\nby hand: '+str(hand_values[1])+' +/- '+str(hand_errors[1]))
    print('\n~~~ Difference in Magnitudes ~~~\nad0101.fits: '+str(differences[0])+'\nad0102.fits: '+str(differences[1]))
    print('\n~~~ Difference in Errors ~~~\nad0101.fits: '+str(error_differences[0])+'\nad0102.fits: '+str(error_differences[1]))
    print('\n~~~ Overall Point Separation ~~~\naphot: '+str(difference_in_points[0])+'\nBy Hand: '+str(difference_in_points[1])+'\n')
    print('\n ~~~ Magnitude Ratios ~~~\nad0101.fits: '+str(mag_ratios[0])+'\nad0102.fits: '+str(mag_ratios[1])+'\n')
    print('\n ~~~ Error Ratios ~~~\nad0101.fits: '+str(err_ratios[0])+'\nad0102.fits: '+str(err_ratios[1])+'\n')
