from imports import *
from data_handling import get_data
from fourier import fourier_fitting

def chi_obs(obj,period,trials):
    times,mags,errors,days=get_data(obj)
    runs=[]
    plt.figure()
    for i in range(trials):
        removed_times,removed_mags,removed_errors=times,mags,errors
        fittable=True
        chis=[fourier_fitting(obj,period,2,5,False,False,20)[2]]
        obs=[len(removed_times)]
        while fittable==True:
            try:
                removed_point=random.randint(0,len(removed_times)-1)
                removed_times=np.delete(removed_times,removed_point)
                removed_mags=np.delete(removed_mags,removed_point)
                removed_errors=np.delete(removed_errors,removed_point)
                



                mean_mag=np.mean(mags)
                mean_mag_error=np.std(mags)/(np.sqrt(len(mags)))
                amp1=1
                period=period
                phi1=np.pi
                amp2=1
                phi2=np.pi

                def fourier_function(t, *params):
                    return params[0]+params[1]*np.sin(2*(np.pi/params[2])*t+params[3])+params[4]*np.sin(5*(np.pi/params[2])*t+params[5])
                    
                fourier_values=[mean_mag,amp1,period,phi1,amp2,phi2]

                amp1_lo=-np.inf
                amp1_hi=np.inf
                amp2_lo=-np.inf
                amp2_hi=np.inf
                p_lo=period/(1+20/100)
                p_hi=period*(1+20/100)
                phi1_lo=0
                phi1_hi=2*np.pi
                phi2_lo=0
                phi2_hi=2*np.pi
                disp_lo=np.mean(mags)-1
                disp_hi=np.mean(mags)+1

                fourier_bounds=([disp_lo,amp1_lo,p_lo,phi1_lo,amp2_lo,phi2_lo],[disp_hi,amp1_hi,p_hi,phi1_hi,amp2_hi,phi2_hi])

                #~~~~~~~~~~~~~ FOURIER FITTING AND PLOT ~~~~~~~~~~~~~

                popt,cov=sp.optimize.curve_fit(fourier_function,removed_times,removed_mags,sigma=removed_errors,p0=fourier_values,bounds=fourier_bounds,check_finite=True,maxfev=10**6)

                #~~~~~~~~~~~~~ CHI SQUARED AND ERRORS FROM CHI SQUARED ~~~~~~~~~~~~~

                def chi_squared(model_params, model, x_data, y_data, y_err):
                    return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
                    
                reduced_chi=chi_squared(popt,fourier_function,times,mags,errors)/len(times)
                obs.append(len(removed_times))
                chis.append(reduced_chi)
            except:
                fittable=False
        runs.append(chis)
        plt.plot(obs,chis,alpha=0.5)
    mean_chis=np.mean(runs,axis=0)
    plt.plot(obs,mean_chis,c='r')
    plt.show()
    
    