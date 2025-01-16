from imports import *
from data_handling import get_data

plt.rcParams.update({'font.size': 22})

def new_fourier(obj,period):
    times,mags,errors,days=get_data(obj)

    possible_harmonics=math.floor(len(times)/2-2)

    if possible_harmonics>6:
        possible_harmonics=6

    print(obj+' Possible Harmonics: '+str(possible_harmonics)+'\n')
    
    def multiharmonics(t, *params):
        period=params[0]
        func=params[1]+params[2]*np.cos((2*np.pi/period)*t+params[3])
        for harmonic in range(2,(int(len(params)/2))):
            func+=params[2*harmonic]*np.cos((2*np.pi/period)*harmonic*t+params[2*harmonic+1])
        return func
    
    def chi_squared(model_params, model, x_data, y_data, y_err):
        return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)

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
        p_lo=period/1.2
        p_hi=period*1.2
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
    
    def folded_multiharmonics(t, *params):
            period=params[0]
            func=params[1]+params[2]*np.cos(2*np.pi*t+params[3])
            for harmonic in range(2,(int(len(params)/2))):
                func+=params[2*harmonic]*np.cos(2*np.pi*harmonic*t+params[2*harmonic+1])
            return func
    
    
    chis=[]
    periods=[]
    harmonic_errors=[]
    asym_errors=[[],[]]
    fitted_params=[]

    foldfig,foldaxs=plt.subplots(3,2)
    #alias_fig,alias_axs=plt.subplots()

    colourtracker=0
    gridtracker=0
    fit_phases=np.linspace(0,1,10000)



    for n in range(1,possible_harmonics+1):
        #fig4, axs4=plt.subplots()
        fit_times=np.linspace(times[0],times[-1],10000)
        #axs4.errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
        harmonic_popt,harmonic_cov=sp.optimize.curve_fit(multiharmonics,times,mags,sigma=errors,p0=generate_harmonic_vals(n),bounds=generate_harmonic_bounds(n),check_finite=True,maxfev=10**6)   
        
        h_degrees_of_freedom=times.size-harmonic_popt.size
        
        periods.append(harmonic_popt[0])
        fitted_period = harmonic_popt[0]
        
        reduced_harmonic_chi=chi_squared(harmonic_popt,multiharmonics,times,mags,errors)/h_degrees_of_freedom
        chis.append(reduced_harmonic_chi)
        
        fitted_params.append(harmonic_popt)

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

        '''plt.figure()
        y,x,_=plt.hist(jackknifed_harmonic_periods,color='b')
        lims=plt.gca().get_xlim()
        x_axis=np.arange(lims[0],lims[1],0.0001)
        plt.axvline(harmonic_negative_percentile,c='r',linestyle='dashed')
        plt.axvline(harmonic_positive_percentile,c='r',linestyle='dashed')
        plt.axvline(fitted_period,c='gold')
        plt.plot(x_axis, (y.max()/100)*norm.pdf(x_axis, fitted_period, harmonic_error_from_jackknifing),c='r',linestyle='dashed') 
        plt.xlabel('Jackknifed Period Fitted by $\chi^2$ Minimisation')
        plt.ylabel('Counts Per Bin')
        plt.show()'''

        harmonic_errors.append(harmonic_error_from_jackknifing)
        asym_errors[0].append(harmonic_negative_error)
        asym_errors[1].append(harmonic_positive_error)

   
        #axs4.plot(fit_times,multiharmonics(fit_times,*harmonic_popt))
        print(f'--------- {n} Harmonic Fit ---------')
        print(f'Period Calculated by {n} Harmonic Fit: {harmonic_popt[0]}')
        print(f'Degrees of Freedom: {h_degrees_of_freedom}')
        print(f'Error Assuming Gaussian Jackknifing Distribution: {harmonic_error_from_jackknifing}')
        print(f'Errors Assuming Non-Gaussian Jackknifing Distribution:\nPositive: {harmonic_positive_error}\nNegative: {harmonic_negative_error}')
        print(f'Reduced Chi Squared for {n} Harmonics: {reduced_harmonic_chi:.2f}\n')


        #axs4.invert_yaxis()
        #axs4.set_ylabel('Magnitude')
        #axs4.set_xlabel('Time (days)')
        #plt.show()

        phases=[]

        for time in times:
            phase=(time%fitted_period)/fitted_period
            phases.append(phase)


        x_pos=gridtracker%2
        y_pos=int((gridtracker-gridtracker%2)/2)

        colours=['r','b','lime','darkviolet','darkorange','gold']
        
        aliasing_periods=np.linspace(period/1.2,period*1.2,1000)
        aliased_chis=[]
        pos=1
        '''
        for aliased_period in aliasing_periods:

            def alias_multiharmonics(t, *params):
                period=aliased_period
                func=params[0]+params[1]*np.cos((2*np.pi/period)*t+params[2])
                for harmonic in range(2,(int(len(params)/2))):
                    func+=params[2*harmonic-1]*np.cos((2*np.pi/period)*harmonic*t+params[2*harmonic])
                return func
            
            def alias_generate_harmonic_vals(n):
                mean_mag=np.mean(mags)
                amp=1
                phi=np.pi

                harmonic_vals=[mean_mag,amp,phi]

                for val in range(2,n+1):
                    harmonic_vals.append(amp)
                    harmonic_vals.append(phi)
                return harmonic_vals
            
            def alias_generate_harmonic_bounds(n):
                amp_lo=-np.inf
                amp_hi=np.inf
                phi_lo=0
                phi_hi=2*np.pi
                disp_lo=np.mean(mags)-1
                disp_hi=np.mean(mags)+1

                harmonic_low_bounds=[disp_lo,amp_lo,phi_lo]
                harmonic_high_bounds=[disp_hi,amp_hi,phi_hi]

                for bound in range(2,n+1):
                    harmonic_low_bounds.append(amp_lo)
                    harmonic_low_bounds.append(phi_lo)
                    harmonic_high_bounds.append(amp_hi)
                    harmonic_high_bounds.append(phi_hi)
                
                return (harmonic_low_bounds,harmonic_high_bounds)
            
            aliased_popt, aliased_cov = sp.optimize.curve_fit(alias_multiharmonics,times,mags,sigma=errors,p0=alias_generate_harmonic_vals(n),bounds=alias_generate_harmonic_bounds(n),check_finite=True,maxfev=10**6)
            
            a_degrees_of_freedom=times.size-aliased_popt.size
            aliased_reduced_chi=chi_squared(aliased_popt,alias_multiharmonics,times,mags,errors)/a_degrees_of_freedom
            delta_chi=aliased_reduced_chi-reduced_harmonic_chi
            print(f'{pos}: {aliased_period} Chi^2={aliased_reduced_chi}')

            aliased_chis.append(delta_chi)

            pos+=1

        alias_axs.plot(aliasing_periods,aliased_chis,c=colours[colourtracker])
        alias_axs.axvline(period/1.2,c='k',linestyle='dashed')
        alias_axs.axvline(period*1.2,c='k',linestyle='dashed')
        alias_axs.axvline(period,c='k',linestyle='dashed')
        alias_axs.axvline(fitted_period,c=colours[colourtracker])
        alias_axs.set_ylabel('$\Delta\chi^2$ of Fourier Decomposition')
        alias_axs.set_xlabel('Period (days)')
        
        '''
        colours=['r','b','lime','darkviolet','darkorange','gold']
        foldaxs[y_pos][x_pos].errorbar(phases,mags,yerr=errors,marker='x',c='k',linestyle='None',capsize=3)
        foldaxs[y_pos][x_pos].plot(fit_phases,folded_multiharmonics(fit_phases,*harmonic_popt),c=colours[colourtracker])
        foldaxs[y_pos][x_pos].invert_yaxis()
        foldaxs[y_pos][x_pos].set_ylabel('Magnitude')
        foldaxs[y_pos][x_pos].set_xlabel('Phase ($\Phi$)')
        foldfig.subplots_adjust(wspace=0, hspace=0)

        if x_pos==1:
            foldaxs[y_pos][x_pos].set_ylabel('')

        colourtracker+=1
        gridtracker+=1
    

    plt.show()



    overfitted=0

    for chi in chis:
        if chi>1:
            overfitted+=1

    fig5,axs5=plt.subplots()
    ns=range(1,possible_harmonics+1)
    axs5.scatter(ns,chis,marker='x',s=60,linestyle='None',c='r')
    axs5.axhline(1,c='b',linestyle='dashed')
    axs5.set_xlabel('Number of fitted Fourier harmonics of the form $m_i\cos{(i\\frac{2\pi}{P}(t-t_0)+\phi_i)}$')
    axs5.set_ylabel('Reduced $\chi^2$ value for fitted data')
    plt.show()

    fig6,axs6=plt.subplots()
    axs6.errorbar(ns,periods,yerr=harmonic_errors,marker='x',c='k',capsize=3,linestyle='None',ecolor='r')
    axs6.errorbar(ns,periods,yerr=asym_errors,marker='x',c='k',capsize=3,linestyle='None',ecolor='b')
    try:
        axs6.axvline(ns[overfitted],c='r',linestyle='dashed')
        axs6.fill_betweenx(axs6.get_ylim(),ns[overfitted],ns[-1]+1,facecolor='r',alpha=.3)
    except:
        pass
    axs6.axhline(period,c='r',linestyle='dashed')
    plt.show()

    max_modes=int(input('\nSelect maximum number of modes: '))

    max_period=periods[max_modes-1]
    max_reduced_chi=chis[max_modes-1]
    max_gaussian_error=harmonic_errors[max_modes-1]
    max_nongaussian_error=[asym_errors[0][max_modes-1],asym_errors[1][max_modes-1]]
    max_params=fitted_params[max_modes-1]

    granularity=10**(-5)

    testing_period=max_period
    delta_chi_1=False

    while delta_chi_1==False:
        testing_period+=granularity

        testing_params=[testing_period]
        other_max_params=max_params[1:]
        for param in other_max_params:
            testing_params.append(param)
        testing_params=np.array(testing_params)
        degrees_of_freedom=times.size-testing_params.size

        testing_chi_squared=chi_squared(testing_params,multiharmonics,times,mags,errors)/degrees_of_freedom

        delta_chi_squared=testing_chi_squared-max_reduced_chi

        if delta_chi_squared>1:
            delta_chi_1=True

    max_positive_chi_1_error=testing_period-max_period

    testing_period=max_period
    delta_chi_1=False

    while delta_chi_1==False:
        testing_period-=granularity

        testing_params=[testing_period]
        other_max_params=max_params[1:]
        for param in other_max_params:
            testing_params.append(param)
        testing_params=np.array(testing_params)
        degrees_of_freedom=times.size-testing_params.size

        testing_chi_squared=chi_squared(testing_params,multiharmonics,times,mags,errors)/degrees_of_freedom

        delta_chi_squared=testing_chi_squared-max_reduced_chi

        if delta_chi_squared>1:
            delta_chi_1=True

    max_negative_chi_1_error=max_period-testing_period
    

    max_phases=[]

    for time in times:
        max_phases.append((time%max_period))

    max_phases=np.array(max_phases)

    def bimodal(t, n1, n2, *params):
        period=params[0]
        f=params[1]+params[2]*np.cos(n1*(np.pi/period)*t+params[3])+params[4]*np.cos(n2*(np.pi/period)*t+params[5])
        return f

    mean_mag=np.mean(mags)
    amp=1
    phi=np.pi

    bimodal_vals=[period,mean_mag,amp,phi,amp,phi]

    amp_lo=-np.inf
    amp_hi=np.inf
    p_lo=period/1.2
    p_hi=period*1.2
    phi_lo=0
    phi_hi=2*np.pi
    disp_lo=np.mean(mags)-1
    disp_hi=np.mean(mags)+1

    bimodal_lo=[p_lo,disp_lo,amp_lo,phi_lo,amp_lo,phi_lo]
    bimodal_hi=[p_hi,disp_hi,amp_hi,phi_hi,amp_hi,phi_hi]

    bimodal_bounds=(bimodal_lo,bimodal_hi)

    n1s=range(1,11)
    n2s=range(1,11)

    colour_plots=np.zeros((10,10))

    for n1 in n1s:
        for n2 in n2s:
            def colour_bimodal(t, *params):
                period=params[0]
                f=params[1]+params[2]*np.cos(n1*(np.pi/period)*t+params[3])+params[4]*np.cos(n2*(np.pi/period)*t+params[5])
                return f

            bimodal_popt,bimodal_cov=sp.optimize.curve_fit(colour_bimodal,max_phases,mags,sigma=errors,p0=bimodal_vals,bounds=bimodal_bounds,check_finite=True,maxfev=10**6)

            bimodal_dof=times.size-bimodal_popt.size


            reduced_bimodal_chi=chi_squared(bimodal_popt,colour_bimodal,max_phases,mags,errors)/bimodal_dof

            colour_plots[n1-1][n2-1]=reduced_bimodal_chi

    plt.figure()
    im=plt.imshow(colour_plots,cmap='plasma',extent=(1,11,1,11),origin='lower',norm='log')
    plt.xlabel('$n_1$ coefficient in 2 mode Fourier function')
    plt.ylabel('$n_2$ coefficient in 2 mode Fourier function')
    plt.colorbar(im,label='Reduced $\chi^2$ value of $n_1$, $n_2$ 2 mode Fourier function')
    plt.axline([1.0,1.0],slope=1,c='k',linestyle='dashed')
    plt.show()

    best_n1=int(input('Best n_1 for bimodal fit: '))
    best_n2=int(input('Best n_2 for bimodal fit: '))

    def best_bimodal(t,*params):
        period=params[0]
        f=params[1]+params[2]*np.cos(best_n1*(np.pi/period)*t+params[3])+params[4]*np.cos(best_n2*(np.pi/period)*t+params[5])
        return f

    best_bimodal_params,best_bimodal_cov=sp.optimize.curve_fit(best_bimodal,times,mags,sigma=errors,p0=bimodal_vals,bounds=bimodal_bounds,check_finite=True,maxfev=10**6)
    best_bimodal_chi=colour_plots[best_n1-1][best_n2-1]
    
    bimodal_period=best_bimodal_params[0]

    jackknifed_bimodal_periods=[]

    for i in range(len(times)):
        jackknifed_mags=np.delete(mags,i)
        jackknifed_times=np.delete(times,i)
        jackknifed_errors=np.delete(errors,i)

        jackknifed_period=sp.optimize.curve_fit(best_bimodal,jackknifed_times,jackknifed_mags,sigma=jackknifed_errors,p0=bimodal_vals,bounds=bimodal_bounds,check_finite=True,maxfev=10**6)[0][0]

        jackknifed_bimodal_periods.append(jackknifed_period)
        
    bimodal_mean_percentile=np.percentile(jackknifed_bimodal_periods,50)
    bimodal_negative_percentile=np.percentile(jackknifed_bimodal_periods,16)
    bimodal_positive_percentile=np.percentile(jackknifed_bimodal_periods,84)

    bimodal_positive_error=bimodal_positive_percentile-bimodal_mean_percentile
    bimodal_negative_error=bimodal_mean_percentile-bimodal_negative_percentile
        
    bimodal_gaussian_error=np.std(jackknifed_bimodal_periods)

    plt.figure()
    plt.hist(jackknifed_bimodal_periods)
    plt.show()

    check_gaussian=bool(input('Gaussian?: '))

    testing_period=bimodal_period
    delta_chi_1=False

    while delta_chi_1==False:
        testing_period+=granularity

        testing_params=[testing_period]
        other_max_params=best_bimodal_params[1:]
        for param in other_max_params:
            testing_params.append(param)
        testing_params=np.array(testing_params)
        degrees_of_freedom=times.size-testing_params.size

        testing_chi_squared=chi_squared(testing_params,best_bimodal,times,mags,errors)/degrees_of_freedom

        delta_chi_squared=testing_chi_squared-best_bimodal_chi

        if delta_chi_squared>1:
            delta_chi_1=True

    bimodal_positive_chi_1_error=testing_period-bimodal_period

    testing_period=bimodal_period
    delta_chi_1=False

    while delta_chi_1==False:
        testing_period-=granularity

        testing_params=[testing_period]
        other_max_params=best_bimodal_params[1:]
        for param in other_max_params:
            testing_params.append(param)
        testing_params=np.array(testing_params)
        degrees_of_freedom=times.size-testing_params.size

        testing_chi_squared=chi_squared(testing_params,best_bimodal,times,mags,errors)/degrees_of_freedom

        delta_chi_squared=testing_chi_squared-best_bimodal_chi

        if delta_chi_squared>1:
            delta_chi_1=True

    bimodal_negative_chi_1_error=bimodal_period-testing_period



    folded_bimodal_params,folded_bimodal_cov=sp.optimize.curve_fit(best_bimodal,max_phases,mags,sigma=errors,p0=bimodal_vals,bounds=bimodal_bounds,check_finite=True,maxfev=10**6)
    
    print(f'~~~ {obj} ~~~\nLiterature Period: {period}\n~~~ HARMONICS ~~~\nMax Number of Modes: {max_modes}\nBest Harmonic Period: {max_period}\nBest Gaussian Error: {max_gaussian_error}\nBest Non-Gaussian Errors :-{max_nongaussian_error[0]},{max_nongaussian_error[1]}\nDelta Chi Squared = 1 Errors: -{max_negative_chi_1_error},{max_positive_chi_1_error}\nHarmonic Reduced Chi Squared: {max_reduced_chi}\n~~~ BIMODAL ~~~\nBimodal Period: {bimodal_period}\nBimodal Gaussian Error {bimodal_gaussian_error}\nBimodal Non-Gaussian Errors :-{bimodal_negative_error},{bimodal_positive_error}\nDelta Chi Squared = 1 Errors: -{bimodal_negative_chi_1_error},{bimodal_positive_chi_1_error}\nBimodal Reduced Chi Squared: {best_bimodal_chi}')

    res_times=np.linspace(0,period,10000)

    res_fig,res_axs=plt.subplots(2,1,height_ratios=(3,1))
    res_axs[0].set_ylabel('Magnitude')
    res_axs[1].set_xlabel('Time (days)')
    res_axs[1].set_ylabel('Norm. Residuals')
    res_axs[0].errorbar(max_phases,mags,yerr=errors,marker='x',c='k',linestyle='None',capsize=3)
    res_axs[0].plot(res_times,multiharmonics(res_times,*max_params),c='r')
    res_axs[0].plot(res_times,bimodal(res_times,best_n1,best_n2,*folded_bimodal_params),c='b')
    res_axs[0].invert_yaxis()
    


    multiharmonic_predictions=[]
    bimodal_predictions=[]
    for phase in max_phases:
        multiharmonic_predictions.append(multiharmonics(phase,*max_params))
        bimodal_predictions.append(bimodal(phase,best_n1,best_n2,*best_bimodal_params))

    multiharmonic_residuals=[]
    bimodal_residuals=[]
    
    pos=0

    for mag in mags:
        multiharmonic_residuals.append((mag-multiharmonic_predictions[pos])/errors[pos])
        bimodal_residuals.append((mag-bimodal_predictions[pos])/errors[pos])
        pos+=1

    res_axs[1].scatter(max_phases,multiharmonic_residuals, c='r',marker='s')
    res_axs[1].scatter(max_phases,bimodal_residuals,c='b',marker='^')
    res_axs[1].axhline(0,c='k',linestyle='dashed')

    res_axs[0].set_xticks([])
    res_axs[1].set_xlim(res_axs[0].get_xlim())
    res_fig.subplots_adjust(hspace=0)

    res_fig.subplots_adjust(right=0.8)
        
    reshist=res_fig.add_axes([0.8,0.11,0.1925,0.1925])

    reshist.hist(bimodal_residuals,color='b',ec='b',orientation='horizontal',alpha=0.5)
    reshist.hist(multiharmonic_residuals,color='r',ec='r',orientation='horizontal',alpha=0.5)
    reshist.set_yticks([])

    plt.savefig('reportfig_'+obj+'.png', dpi = 300,bbox_inches='tight')
    plt.show()


    if check_gaussian==True:
        fit_errors=[]
        if bimodal_gaussian_error>abs(bimodal_negative_chi_1_error):
            fit_errors.append(bimodal_gaussian_error)
        else:
            fit_errors.append(bimodal_negative_chi_1_error)

        if bimodal_gaussian_error>abs(bimodal_positive_chi_1_error):
            fit_errors.append(bimodal_gaussian_error)
        else:
            fit_errors.append(bimodal_positive_chi_1_error)
    else:
        fit_errors=[]
        if bimodal_negative_error>abs(bimodal_negative_chi_1_error):
            fit_errors.append(bimodal_negative_error)
        else:
            fit_errors.append(bimodal_negative_chi_1_error)

        if bimodal_positive_error>abs(bimodal_positive_chi_1_error):
            fit_errors.append(bimodal_positive_error)
        else:
            fit_errors.append(bimodal_positive_chi_1_error)
    
    m_0=np.median(mags)
    m_0_err=np.std(mags)/(np.sqrt(len(mags)))   

    fitting_params=[obj, period, max_modes, max_period, max_reduced_chi, best_n1, best_n2, bimodal_period, fit_errors, best_bimodal_chi]
    
    print(f'\n{obj}\nLiterature Period: {period}\n')
    print(f'Maximum Number of Fourier Modes Used in Full Decomposition: {max_modes}\n')
    print(f'Fourier Decomposition Period: {max_period}\n')
    print(f'Errors from chisq: -{max_negative_chi_1_error},+{max_positive_chi_1_error}\n')
    print(f'Errors from Jackkinfing Assuming Gaussian Distribution: {max_gaussian_error}\n')
    print(f'Errors from Jackknifing Assuming Non-Gaussian Distribution: -{max_nongaussian_error[0]},+{max_nongaussian_error[1]}\n')
    print(f'Reduced chisq of Fourier Decomposition: {max_reduced_chi}\n')

    print(f'Fourier Coefficients Used in Bimodal Model: {best_n1}, {best_n2}\n')
    print(f'Bimodal Period: {bimodal_period}\n')
    print(f'Errors from chisq: -{bimodal_negative_chi_1_error},+{bimodal_positive_chi_1_error}\n')
    print(f'Errors from Jackkinfing Assuming Gaussian Distribution: {bimodal_gaussian_error}\n')
    print(f'Errors from Jackknifing Assuming Non-Gaussian Distribution: -{bimodal_negative_error},+{bimodal_positive_error}\n')
    print(f'Reduced chisq of Fourier Decomposition: {best_bimodal_chi}\n')

    print(f'delta chisq={best_bimodal_chi-max_reduced_chi}\n')

    print(f'Average Magnitude: {m_0}pm{m_0_err}\n')

    plot_data=''

    plot_data+=f'\n{obj}\nLiterature Period: {period}\n'
    plot_data+=f'Maximum Number of Fourier Modes Used in Full Decomposition: {max_modes}\n'
    plot_data+=f'Fourier Decomposition Period: {max_period}\n'
    plot_data+=f'Errors from chisq: -{max_negative_chi_1_error},+{max_positive_chi_1_error}\n'
    plot_data+=f'Errors from Jackkinfing Assuming Gaussian Distribution: {max_gaussian_error}\n'
    plot_data+=f'Errors from Jackknifing Assuming Non-Gaussian Distribution: -{max_nongaussian_error[0]},+{max_nongaussian_error[1]}\n'
    plot_data+=f'Reduced chisq of Fourier Decomposition: {max_reduced_chi}\n'

    plot_data+=f'Fourier Coefficients Used in Bimodal Model: {best_n1}, {best_n2}\n'
    plot_data+=f'Bimodal Period: {bimodal_period}\n'
    plot_data+=f'Errors from chisq: -{bimodal_negative_chi_1_error},+{bimodal_positive_chi_1_error}\n'
    plot_data+=f'Errors from Jackkinfing Assuming Gaussian Distribution: {bimodal_gaussian_error}\n'
    plot_data+=f'Errors from Jackknifing Assuming Non-Gaussian Distribution: -{bimodal_negative_error},+{bimodal_positive_error}\n'
    plot_data+=f'Reduced chisq of Fourier Decomposition: {best_bimodal_chi}\n'

    plot_data+=f'delta chisq={best_bimodal_chi-max_reduced_chi}\n'

    plot_data+=f'Average Magnitude: {m_0}pm{m_0_err}\n'


    f=open('plot_data.txt','a')
    f.write(plot_data)
    f.close()

    return m_0,m_0_err,bimodal_period,fit_errors,fitting_params
    
    


        
    




    
    

    

    



