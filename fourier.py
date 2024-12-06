from imports import *
from data_handling import get_data


def fourier_fitting(obj,period,n1,n2,show_plots,folded_dates,bound_percentage):
    #~~~~~~~~~~~~~ DATA IMPORT ~~~~~~~~~~~~~
    times,mags,errors,days=get_data(obj)

    #~~~~~~~~~~~~~ FITTING MODEL, PARAMETERS AND BOUNDS ~~~~~~~~~~~~~
    mean_mag=np.mean(mags)
    mean_mag_error=np.std(mags)/(np.sqrt(len(mags)))
    amp1=1
    period=period
    phi1=np.pi
    amp2=1
    phi2=np.pi

    def fourier_function(t, *params):
        return params[0]+params[1]*np.cos(n1*(np.pi/params[2])*t+params[3])+params[4]*np.cos(n2*(np.pi/params[2])*t+params[5])
    
    fourier_values=[mean_mag,amp1,period,phi1,amp2,phi2]

    amp1_lo=-np.inf
    amp1_hi=np.inf
    amp2_lo=-np.inf
    amp2_hi=np.inf
    p_lo=period/(1+bound_percentage/100)
    p_hi=period*(1+bound_percentage/100)
    phi1_lo=0
    phi1_hi=2*np.pi
    phi2_lo=0
    phi2_hi=2*np.pi
    disp_lo=np.mean(mags)-1
    disp_hi=np.mean(mags)+1

    fourier_bounds=([disp_lo,amp1_lo,p_lo,phi1_lo,amp2_lo,phi2_lo],[disp_hi,amp1_hi,p_hi,phi1_hi,amp2_hi,phi2_hi])

    #~~~~~~~~~~~~~ FOURIER FITTING ~~~~~~~~~~~~~

    popt,cov=sp.optimize.curve_fit(fourier_function,times,mags,sigma=errors,p0=fourier_values,bounds=fourier_bounds,check_finite=True,maxfev=10**6)

    smooth_x=np.linspace(times[0], times[-1], 1000)

    fitted_period=popt[2]

    #~~~~~~~~~~~~~ RESIDUALS ~~~~~~~~~~~~~

    predictions=[]
    for i in times:
        predictions.append(fourier_function(i,*popt))
    
    residuals=[]
    pos=0
    for val in mags:
        residuals.append((val-predictions[pos])/errors[pos])
        pos+=1
    
    #~~~~~~~~~~~~~ INITIAL PLOTTING ~~~~~~~~~~~~~
    output_popt=np.round(popt,2)
    if show_plots==True:
        fig1, axs1 = plt.subplots(2,1,height_ratios=(3,1))
        axs1[0].errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
        axs1[0].plot(smooth_x,fourier_function(smooth_x, *popt),c='r')
        axs1[0].set_xlabel('Time (days)') 
        axs1[0].set_ylabel('Magnitude')
        axs1[0].invert_yaxis()
        axs1[0].set_xticks([])
        axs1[1].scatter(times,residuals,c='k',marker='x')
        axs1[1].set_xlabel('Time (days)')
        axs1[1].set_ylabel('Normalised Residuals')

        limit=max(abs(np.array(residuals)))+1
        axs1[1].set_ylim([-limit,limit])
        #axs1[1].axhline(0,c='r',linestyle='dashed',alpha=0.6)
        ticks1=axs1[1].get_yticks()/2
        pos=0
        for i in ticks1:
            if pos%2==1:
                axs1[1].axhline(i,c='k',linestyle='dashed',alpha=0.6)
            else:
                axs1[1].axhline(i,c='r',linestyle='dashed',alpha=0.6)
            pos+=1

        fig1.subplots_adjust(right=0.8)
        
        fig1hist=fig1.add_axes([0.8,0.11,0.1925,0.1925])

        fig1hist.hist(residuals, orientation='horizontal')

        fig1hist.set_ylim([-limit,limit])
        fig1hist.set_yticks([])

        fig1.subplots_adjust(hspace=0)
        plt.show()

    #~~~~~~~~~~~~~ JACKKNIFING ~~~~~~~~~~~~~

    jackknifed_periods=[]

    for i in range(len(times)):
        jackknifed_mags=np.delete(mags,i)
        jackknifed_times=np.delete(times,i)
        jackknifed_errors=np.delete(errors,i)

        jackknifed_period=sp.optimize.curve_fit(fourier_function,jackknifed_times,jackknifed_mags,sigma=jackknifed_errors,p0=fourier_values,bounds=fourier_bounds,check_finite=True,maxfev=10**6)[0][2]

        jackknifed_periods.append(jackknifed_period)

    if show_plots==True:
        plt.figure()
        plt.hist(jackknifed_periods,bins=100)
        plt.xlabel('Jackknifed Period (days)')
        plt.show()
        

    mean_percentile=np.percentile(jackknifed_periods,50)
    negative_percentile=np.percentile(jackknifed_periods,16)
    positive_percentile=np.percentile(jackknifed_periods,84)

    positive_error=positive_percentile-mean_percentile
    negative_error=mean_percentile-negative_percentile
    
    error_from_jackknifing=np.std(jackknifed_periods)

    #~~~~~~~~~~~~~ CHI SQUARED AND ERRORS FROM CHI SQUARED ~~~~~~~~~~~~~

    def chi_squared(model_params, model, x_data, y_data, y_err):
        return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
    degrees_of_freedom=times.size-popt.size
    reduced_chi=chi_squared(popt,fourier_function,times,mags,errors)/degrees_of_freedom

    min_period=period*(1+bound_percentage/100)
    max_period=period/(1+bound_percentage/100)

    periods=np.linspace(min_period,max_period,10000)
    fourier_chis=[]

    for a_period in periods:

        a_mean_mag=np.mean(mags)
        a_amp1=1
        a_period=a_period
        a_phi1=np.pi
        a_amp2=1
        a_phi2=np.pi

        def alias_fourier_function(t, *params):
            return params[0]+params[1]*np.cos(2*(np.pi/a_period)*t+params[2])+params[3]*np.cos(5*(np.pi/period)*t+params[4])

        a_fourier_values=[a_mean_mag,a_amp1,a_phi1,a_amp2,a_phi2]

        a_fourier_popt,a_fourier_cov=sp.optimize.curve_fit(alias_fourier_function,times,mags,sigma=errors,p0=a_fourier_values,check_finite=True,maxfev=10**6)
    

        a_fourier_chi_val=chi_squared(a_fourier_popt,alias_fourier_function,times,mags,errors)/degrees_of_freedom

        fourier_chis=np.append(fourier_chis,a_fourier_chi_val)
    
    fourier_bound_period=popt[2]
    bound_error=error_from_jackknifing
    delta_chis=fourier_chis-min(fourier_chis)

    min_index=np.argmin(delta_chis)

    low_pos=min_index
    found_low_err=False
    while found_low_err==False:
        if delta_chis[low_pos]>1:
              found_low_err=True
        else:
              low_pos-=1
    
    high_pos=min_index
    found_high_err=False
    while found_high_err==False:
        if delta_chis[high_pos]>1:
              found_high_err=True
        else:
              high_pos+=1
    
    low_err=periods[min_index]-periods[low_pos]
    high_err=periods[high_pos]-periods[min_index]

    if show_plots==True:
        fig,ax=plt.subplots()
        plt.plot(periods,delta_chis,c='r')
        ax.axvline(fourier_bound_period,c='b',linestyle='dashed')
        ax.axvline(fourier_bound_period-bound_error,c='b',linestyle='dashed',alpha=0.5)
        ax.axvline(fourier_bound_period+bound_error,c='b',linestyle='dashed',alpha=0.5)
        ax.axhline(1,c='g',linestyle='dashed')
        ax.axvline(min_period,c='k',linestyle='dashed')
        ax.axvline(max_period,c='k',linestyle='dashed')
        ax.axvline(fourier_bound_period-low_err,c='g',linestyle='dashed',alpha=.5)
        ax.axvline(fourier_bound_period+high_err,c='g',linestyle='dashed',alpha=.5)
        plt.ylabel('$\Delta\chi^2$')
        plt.xlabel('Period (days)')
        plt.show()

    #~~~~~~~~~~~~~ DATA READOUTS ~~~~~~~~~~~~~

    if show_plots==True:
        print('Period (days): '+str(popt[2]))
        print('Positive Error from Jackknifing: '+str(positive_error))
        print('Bottom Error from Jackknifing: '+str(negative_error))
        print('Error from Standard Deviation of Jackknifing: '+str(error_from_jackknifing))
        print('Error from Chi Squared + 1: '+str(np.mean([abs(low_err),abs(high_err)])))
        print('Fitted Function: '+str(output_popt[0])+'+'+str(output_popt[1])+'cos(2*pi/'+str(output_popt[2])+'t+'+str(output_popt[3])+')+'+str(output_popt[4])+'cos(5*pi/'+str(output_popt[2])+'t+'+str(output_popt[5])+')')
        print('Degrees of Freedom: '+str(degrees_of_freedom))
        print('Reduced Chi Squared: '+str(reduced_chi))

        error_up=[popt[0],popt[1],popt[2]+error_from_jackknifing,popt[3],popt[4],popt[5]]

        error_down=[popt[0],popt[1],popt[2]-error_from_jackknifing,popt[3],popt[4],popt[5]]

        percent_error_up=[popt[0],popt[1],popt[2]+positive_error,popt[3],popt[4],popt[5]]

        percent_error_down=[popt[0],popt[1],popt[2]-negative_error,popt[3],popt[4],popt[5]]

        chi_error_up=[popt[0],popt[1],popt[2]+high_err,popt[3],popt[4],popt[5]]

        chi_error_down=[popt[0],popt[1],popt[2]-low_err,popt[3],popt[4],popt[5]]

        fig2,axs2=plt.subplots(2,1,height_ratios=(3,1))
        axs2[0].errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *popt),c='r')
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *error_up),c='r',linestyle='dashed',alpha=0.5)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *error_down),c='r',linestyle='dashed',alpha=0.5)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *percent_error_up),c='b',linestyle='dashed',alpha=0.5)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *percent_error_down),c='b',linestyle='dashed',alpha=0.5)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *chi_error_up),c='g',linestyle='dashed',alpha=0.5)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *chi_error_down),c='g',linestyle='dashed',alpha=0.5)
        axs2[0].set_xlabel('Time (days)') 
        axs2[0].set_ylabel('Magnitude')
        axs2[0].invert_yaxis()
        axs2[0].set_xticks([])

        axs2[1].scatter(times,residuals,c='k',marker='x')
        axs2[1].set_xlabel('Time (days)')
        axs2[1].set_ylabel('Normalised Residuals')

        limit=max(abs(np.array(residuals)))+1
        axs2[1].set_ylim([-limit,limit])
        ticks2=axs2[1].get_yticks()/2
        pos=0
        for i in ticks2:
            if pos%2==1:
                axs2[1].axhline(i,c='k',linestyle='dashed',alpha=0.6)
            else:
                axs2[1].axhline(i,c='r',linestyle='dashed',alpha=0.6)
            pos+=1

        fig2.subplots_adjust(right=0.8)
        
        fig2hist=fig2.add_axes([0.8,0.11,0.1925,0.1925])

        fig2hist.hist(residuals, orientation='horizontal')

        fig2hist.set_ylim([-limit,limit])
        fig2hist.set_yticks([])

        

        fig2.subplots_adjust(hspace=0)
        plt.show()


    #~~~~~~~~~~~~~ FOLDING ~~~~~~~~~~~~~

    folded_times=[]
    for i in times:
        folded_times.append(i)
    folded_phase=[]
    for i in range(len(folded_times)):
        folded_times[i]=folded_times[i]%fitted_period
        folded_phase.append(folded_times[i]/fitted_period)

    folded_fit_times=np.linspace(0,fitted_period,1000)
    
    days=np.array(days)

    mjd_dates=astro.time.Time(days,format='mjd')

    utc=mjd_dates.to_datetime()

    currentdate=astro.time.Time.now()
    currentdate=astro.time.Time(currentdate)
    currentdate-=times[0]
    folded_currentdate=currentdate.mjd%fitted_period

    folded_predictions=[]
    for i in folded_times:
        folded_predictions.append(fourier_function(i,*popt))
    
    folded_residuals=[]
    pos=0
    for val in mags:
        folded_residuals.append((val-folded_predictions[pos])/errors[pos])
        pos+=1
    
    date_labels=[]

    for i in utc:
        date_labels.append(str(i.day)+'/'+str(i.month)+'/'+str(i.year))    
    
    if show_plots==True:

        fig3,axs3=plt.subplots(2,1,height_ratios=(3,1))
        axs3[0].errorbar(folded_times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
        axs3[0].plot(folded_fit_times,fourier_function(folded_fit_times,*popt),c='r')
        axs3[0].plot(folded_fit_times,fourier_function(folded_fit_times,*error_up),c='r',linestyle='dashed',alpha=0.5)
        axs3[0].plot(folded_fit_times,fourier_function(folded_fit_times,*error_down),c='r',linestyle='dashed',alpha=0.5)
        axs3[0].plot(folded_fit_times,fourier_function(folded_fit_times, *percent_error_up),c='b',linestyle='dashed',alpha=0.5)
        axs3[0].plot(folded_fit_times,fourier_function(folded_fit_times, *percent_error_down),c='b',linestyle='dashed',alpha=0.5)
        axs3[0].plot(folded_fit_times,fourier_function(folded_fit_times, *chi_error_up),c='g',linestyle='dashed',alpha=0.5)
        axs3[0].plot(folded_fit_times,fourier_function(folded_fit_times, *chi_error_down),c='g',linestyle='dashed',alpha=0.5)
        axs3[0].set_xticks([])
        if folded_dates==True:
            for i in range(len(date_labels)):
                axs3[0].text(folded_times[i],mags[i],date_labels[i])
            axs3[0].axvline(folded_currentdate,c='b',linestyle='dashed')
        axs3[0].set_ylabel('Magnitude')
        axs3[0].set_xlabel('Time in Period (days)')
        axs3[0].invert_yaxis()


        axs3[1].scatter(folded_times,folded_residuals,c='k',marker='x')
        axs3[1].set_xlabel('Time (days)')
        axs3[1].set_ylabel('Normalised Residuals')

        limit=max(abs(np.array(folded_residuals)))+1
        axs3[1].set_ylim([-limit,limit])

        axs3[1].set_ylim([-limit,limit])
        ticks3=axs3[1].get_yticks()/2
        pos=0
        for i in ticks3:
            if pos%2==1:
                axs3[1].axhline(i,c='k',linestyle='dashed',alpha=0.6)
            else:
                axs3[1].axhline(i,c='r',linestyle='dashed',alpha=0.6)
            pos+=1

        fig3.subplots_adjust(right=0.8)
        
        fig3hist=fig3.add_axes([0.8,0.11,0.1925,0.1925])

        fig3hist.hist(folded_residuals, orientation='horizontal')

        fig3hist.set_ylim([-limit,limit])
        fig3hist.set_yticks([])

        fig3.subplots_adjust(hspace=0)

        plt.show()

    #~~~~~~~~~~~~~ MULTIPLE HARMONICS ~~~~~~~~~~~~~

    possible_harmonics=math.floor(len(times)/2-2)
    print(obj+' Possible Harmonics: '+str(possible_harmonics))
    
    def multiharmonics(t, *params):
        period=params[0]
        func=params[1]+params[2]*np.cos((2*np.pi/period)*t+params[3])
        for harmonic in range(2,(int(len(params)/2))):
            func+=params[2*harmonic]*np.cos((2*np.pi/period)*harmonic*t+params[2*harmonic+1])
        return func

    def generate_harmonic_vals(n):
        mean_mag=np.mean(mags)
        amp=1
        phi=np.pi

        harmonic_vals=[period,mean_mag,amp,phi]

        for val in range(2,n+1):
            harmonic_vals.append(amp)
            harmonic_vals.append(phi)
        return harmonic_vals
    
    def generate_harmonic_bounds(n):
        amp_lo=-np.inf
        amp_hi=np.inf
        p_lo=period/(1+bound_percentage/100)
        p_hi=period*(1+bound_percentage/100)
        phi_lo=0
        phi_hi=2*np.pi
        disp_lo=np.mean(mags)-1
        disp_hi=np.mean(mags)+1

        harmonic_low_bounds=[p_lo,disp_lo,amp_lo,phi_lo]
        harmonic_high_bounds=[p_hi,disp_hi,amp_hi,phi_hi]

        for bound in range(2,n+1):
            harmonic_low_bounds.append(amp_lo)
            harmonic_low_bounds.append(phi_lo)
            harmonic_high_bounds.append(amp_hi)
            harmonic_high_bounds.append(phi_hi)
        
        return (harmonic_low_bounds,harmonic_high_bounds)
    

    if show_plots==True:
        fig4, axs4=plt.subplots(2,1)
        fit_times=np.linspace(times[0],times[-1],1000)
        axs4[0].errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
        axs4[0].plot(fit_times,fourier_function(fit_times,*popt),c='r')
        axs4[1].errorbar(folded_times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
        axs4[1].plot(folded_fit_times,fourier_function(folded_fit_times,*popt),c='r')

    chis=[reduced_chi]
    harmonic_periods=[fitted_period]
    harmonic_period_errors=[error_from_jackknifing]
    harmonics=range(2,possible_harmonics+1)

    for n in range(3,possible_harmonics+1):
        harmonic_popt,harmonic_cov=sp.optimize.curve_fit(multiharmonics,times,mags,sigma=errors,p0=generate_harmonic_vals(n),bounds=generate_harmonic_bounds(n),check_finite=True,maxfev=10**6)   
        reduced_harmonic_chi=chi_squared(harmonic_popt,multiharmonics,times,mags,errors)/degrees_of_freedom

        chis.append(reduced_harmonic_chi)
        harmonic_periods.append(harmonic_popt[0])

        jackknifed_harmonic_periods=[]

        for i in range(len(times)):
            jackknifed_mags=np.delete(mags,i)
            jackknifed_times=np.delete(times,i)
            jackknifed_errors=np.delete(errors,i)

            jackknifed_period=sp.optimize.curve_fit(multiharmonics,jackknifed_times,jackknifed_mags,sigma=jackknifed_errors,p0=generate_harmonic_vals(n),bounds=generate_harmonic_bounds(n),check_finite=True,maxfev=10**6)[0][0]

            jackknifed_harmonic_periods.append(jackknifed_period)
        
        harmonic_mean_percentile=np.percentile(jackknifed_harmonic_periods,50)
        harmonic_negative_percentile=np.percentile(jackknifed_harmonic_periods,16)
        harmonic_positive_percentile=np.percentile(jackknifed_harmonic_periods,84)

        harmonic_positive_error=harmonic_positive_percentile-harmonic_mean_percentile
        harmonic_negative_error=harmonic_mean_percentile-harmonic_negative_percentile
        
        harmonic_error_from_jackknifing=np.std(jackknifed_harmonic_periods)



        if show_plots==True:
            if reduced_harmonic_chi>1:
                axs4[0].plot(fit_times,multiharmonics(fit_times,*harmonic_popt),alpha=0.5)
                print(f'Period Calculated by {n} Harmonic Fit: {harmonic_popt[0]}')
                print(f'Error Assuming Gaussian Jackknifing Distribution: {harmonic_error_from_jackknifing}')
                print(f'Errors Assuming Non-Gaussian Jackknifing Distribution:\nPositive: {harmonic_positive_error}\nNegative: {harmonic_negative_error}')
                print(f'Reduced Chi Squared for {n} Harmonics: {reduced_harmonic_chi:.2f}')
                axs4[1].plot(folded_fit_times,multiharmonics(folded_fit_times,*harmonic_popt),label='n='+str(n),alpha=0.5)

    if show_plots==True:
        axs4[1].legend()
        axs4[0].invert_yaxis()
        axs4[1].invert_yaxis()
        plt.show()

        fig5,axs5=plt.subplots()
        ns=range(2,possible_harmonics+1)
        axs5.scatter(ns,chis,marker='x',linestyle='None',c='k')
        axs5.axhline(1,c='r',linestyle='dashed')
        axs5.set_xlabel('Number of fitted Fourier harmonics of the form $A_n\cos(n\\frac{2\pi}{T}t+\phi_n)$')
        axs5.set_ylabel('Reduced $\chi^2$ value for fitted data')
        plt.show()


    return fitted_period,error_from_jackknifing,reduced_chi,popt[0], mean_mag_error
