from imports import *
from data_handling import get_data

def modes_plot(obj,period,modes,show_plots):
    times,mags,errors,days=get_data(obj)
    mode_vals=np.linspace(1,modes,modes)
    chis=[]
    for totalmode in range(1,modes+1):
            mean_mag=np.mean(mags)
            amp1=1
            period=period
            phi1=np.pi
                
            initial_values=[mean_mag,amp1,period,phi1]
            for i in range(2,totalmode+1):
                initial_values.append(1)
                initial_values.append(np.pi)

            def fourier_function(t,*params):
                f = params[0]+params[1]*np.sin(1*(np.pi/params[2])*t+params[3])
                for mode in range(2,totalmode+1):
                     f+=params[2*mode]*np.sin(mode*(np.pi/params[2])*t+params[2*mode+1])
                return f
            
            amp_lo=-np.inf
            p_lo=period/(1+20/100)
            phi_lo=0
            disp_lo=np.mean(mags)-1

            lower_bound=[disp_lo,amp_lo,p_lo,phi_lo]

            for i in range(2,totalmode+1):
                 lower_bound.append(amp_lo)
                 lower_bound.append(phi_lo)

            amp_hi=np.inf
            p_hi=period*(1+20/100)
            phi_hi=np.pi
            disp_hi=np.mean(mags)+1

            upper_bound=[disp_hi,amp_hi,p_hi,phi_hi]

            for i in range(2,totalmode+1):
                 upper_bound.append(amp_hi)
                 upper_bound.append(phi_hi)
            
            bounds=(lower_bound,upper_bound)

            popt,cov=sp.optimize.curve_fit(fourier_function,times,mags,sigma=errors,p0=initial_values,bounds=bounds,check_finite=True,maxfev=10**6)

            def chi_squared(model_params, model, x_data, y_data, y_err):
                return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
            degrees_of_freedom=times.size-popt.size
            reduced_chi=chi_squared(popt,fourier_function,times,mags,errors)/degrees_of_freedom
            chis.append(reduced_chi)

    overfitted=False
    pos=0
    while overfitted!=True:
        if chis[pos]<1:
            overfitted=True
        else:
            pos+=1
    overfitted_mode=mode_vals[pos]
    if show_plots==True:
        plt.figure()
        plt.plot(mode_vals,chis,c='r')
        plt.xlabel('Number of fitted fourier modes of the form $A_n\sin(n\\frac{\pi}{T}t+\phi_n)$')
        plt.ylabel('Reduced $\chi^2$ value for fitted data')
        plt.axvline(mode_vals[pos],c='b',linestyle='dashed')
        plt.show()
    return overfitted_mode

                