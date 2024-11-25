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
        return params[0]+params[1]*np.sin(n1*(np.pi/params[2])*t+params[3])+params[4]*np.sin(n2*(np.pi/params[2])*t+params[5])
    
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


        fig1.subplots_adjust(hspace=0)
        plt.show()

    #~~~~~~~~~~~~~ CHI SQUARED AND ERRORS FROM CHI SQUARED ~~~~~~~~~~~~~

    def chi_squared(model_params, model, x_data, y_data, y_err):
        return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
    degrees_of_freedom=times.size-popt.size
    reduced_chi=chi_squared(popt,fourier_function,times,mags,errors)/degrees_of_freedom

    granularity=0.001
    error_period=popt[2]
    current_chi=reduced_chi

    while current_chi<reduced_chi+1:
        error_period+=granularity
        current_chi=chi_squared([popt[0],popt[1],error_period,popt[3],popt[4],popt[5]],fourier_function,times,mags,errors)

    upper_error=abs(popt[2]-error_period)

    error_period=popt[2]
    current_chi=reduced_chi

    while current_chi<reduced_chi+1:
        error_period-=granularity
        current_chi=chi_squared([popt[0],popt[1],error_period,popt[3],popt[4],popt[5]],fourier_function,times,mags,errors)

    lower_error=abs(popt[2]-error_period)

    fitted_period=popt[2]
    chi_plus_1_error=np.mean([lower_error,upper_error])

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

    #~~~~~~~~~~~~~ DATA READOUTS ~~~~~~~~~~~~~

    if show_plots==True:
        print('Period (days): '+str(popt[2]))
        #print('Error From Chi+1: '+str(chi_plus_1_error))
        print('Positive Error from Jackknifing: '+str(positive_error))
        print('Bottom Error from Jackknifing: '+str(negative_error))
        print('Error from Standard Deviation of Jackknifing: '+str(error_from_jackknifing))
        print('Fitted Function: '+str(output_popt[0])+'+'+str(output_popt[1])+'sin(2*pi/'+str(output_popt[2])+'t+'+str(output_popt[3])+')+'+str(output_popt[4])+'sin(5*pi/'+str(output_popt[2])+'t+'+str(output_popt[5])+')')
        print('Degrees of Freedom: '+str(degrees_of_freedom))
        print('Reduced Chi Squared: '+str(reduced_chi))

        error_up=[popt[0],popt[1],popt[2]+error_from_jackknifing,popt[3],popt[4],popt[5]]

        error_down=[popt[0],popt[1],popt[2]-error_from_jackknifing,popt[3],popt[4],popt[5]]

        percent_error_up=[popt[0],popt[1],popt[2]+positive_error,popt[3],popt[4],popt[5]]

        percent_error_down=[popt[0],popt[1],popt[2]-negative_error,popt[3],popt[4],popt[5]]

        fig2,axs2=plt.subplots(2,1,height_ratios=(3,1))
        axs2[0].errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *popt),c='r')
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *error_up),c='r',linestyle='dashed',alpha=0.5)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *error_down),c='r',linestyle='dashed',alpha=0.5)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *percent_error_up),c='b',linestyle='dashed',alpha=0.5)
        axs2[0].plot(smooth_x,fourier_function(smooth_x, *percent_error_down),c='b',linestyle='dashed',alpha=0.5)
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
        fig2.subplots_adjust(hspace=0)
        plt.show()


    #~~~~~~~~~~~~~ FOLDING ~~~~~~~~~~~~~

    folded_times=times
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

        fig3.subplots_adjust(hspace=0)

        plt.show()

    return fitted_period,error_from_jackknifing,reduced_chi,popt[0], mean_mag_error
