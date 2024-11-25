from imports import *
plt.rcParams.update({'font.size': 22})


def get_periods():
    
    df = pd.read_csv("periodmagdata.txt", delimiter=" ", header=None)
    
    mags = df[0].to_numpy()
    mags_errs = df[1].to_numpy()
    periods = df[2].to_numpy()
    periods_errs = df[3].to_numpy()
    
    
    return mags, mags_errs, periods, periods_errs


def get_parallaxes():
    
    objs=next(os.walk('.'))[1]
    objs = objs[1:-1]
    
    gaiadr3_ids = []
    
    Simbad.add_votable_fields('ids')
    
    for obj in objs:
        
        result = Simbad.query_object(obj)
        
        for id in result['ids'][0].split('|'):
            
            if 'Gaia DR3' in id:
                
                gaia_id = id.split(' ')[-1]
            
                break
            
        gaiadr3_ids = np.append(gaiadr3_ids,gaia_id)
            
    
    
    
        
    """query = f"SELECT parallax, main_id \
    FROM gaiadr2.gaia_source \
    WHERE source_id = " + gaiadr3_ids[0] + "\
    OR source_id = " + gaiadr3_ids[1] + "\
    OR source_id = " + gaiadr3_ids[2] + "\
    OR source_id = " + gaiadr3_ids[3] + "\
    OR source_id = " + gaiadr3_ids[4] + "\
    OR source_id = " + gaiadr3_ids[5] + "\
    OR source_id = " + gaiadr3_ids[6] + "\
    OR source_id = " + gaiadr3_ids[7] + "\
    OR source_id = " + gaiadr3_ids[8] + "\
    OR source_id = " + gaiadr3_ids[9] + "\
    OR source_id = " + gaiadr3_ids[10] + "\
    OR source_id = " + gaiadr3_ids[11] + "\
    OR source_id = " + gaiadr3_ids[12] + "\
    OR source_id = " + gaiadr3_ids[13] + "\
    OR source_id = " + gaiadr3_ids[14] + "\
    OR source_id = " + gaiadr3_ids[15] + "\
    OR source_id = " + gaiadr3_ids[16] + "\
    OR source_id = " + gaiadr3_ids[17] + "\
    OR source_id = " + gaiadr3_ids[18] + "\
    OR source_id = " + gaiadr3_ids[19] + "\
    OR source_id = " + gaiadr3_ids[20] + "\
    OR source_id = " + gaiadr3_ids[21] + "\
    OR source_id = " + gaiadr3_ids[22] + "\
    OR source_id = " + gaiadr3_ids[23]"""

    parallaxes = []
    parallaxes_errs = []

    for i in range(len(objs)):
        
        query = f"SELECT parallax, parallax_error \
        FROM gaiadr3.gaia_source \
        WHERE source_id = " + gaiadr3_ids[i]

        job     = Gaia.launch_job_async(query)
        results = job.get_results()
        
        parallaxes = np.append(parallaxes,np.float64(results['parallax']))
        parallaxes_errs = np.append(parallaxes_errs,np.float64(results['parallax_error']))
        
    dist_pc = 1/(0.001 * parallaxes)                                        #GAIA parallaxes are in mas, convert to arcseconds here7
    
    dist_pc_errs = (0.001*parallaxes_errs)/((0.001*parallaxes)**2)
    
    df = pd.DataFrame({"Distance" : dist_pc,"Error" : dist_pc_errs})
    
    df.to_csv('parallax.csv',sep=' ')
    
    return dist_pc, dist_pc_errs


def read_parallaxes():
    
    df = pd.read_csv("parallax.csv",delimiter=" ")
    
    dist_pc = df['Distance'].to_numpy()
    dist_pc_errs = df['Error'].to_numpy()
    
    df = pd.read_csv("ddodata.txt", header=None, delimiter=" ")
    
    bvcorrection = df[0].to_numpy()
    
    return dist_pc, dist_pc_errs, bvcorrection


def fitting_model(x, params):
    return params[0]*np.log10(x) + params[1]


def chi_squared(model_params, model, x_data, y_data, y_err):
    chi_sq = np.sum(((y_data - model(x_data, model_params))/y_err)**2)
    print(chi_sq)
    return chi_sq

def mean_absolute_deviation(model_params, model, x_data, y_data, npoints):
    mad = np.sum(np.abs(y_data - model(x_data, model_params)))/npoints
    print(mad)
    return mad


def plot_calc(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection):
    
    objs=next(os.walk('.'))[1]
    objs = np.array(objs[1:-1])
    
    print(bvcorrection)
    
    abs_mags = (mags - 3.1* bvcorrection) - 5*np.log10(dist_pc) + 5
    
    abs_mags_errs = np.sqrt(mags_errs**2 + (5*(dist_pc_errs/(np.log(10)*dist_pc)))**2)
    
    sort_index = np.argsort(periods)
    
    abs_mags = abs_mags[sort_index]
    periods = periods[sort_index]
    abs_mags_errs = abs_mags_errs[sort_index]
    periods_errs = periods_errs[sort_index]
    objs = objs[sort_index]
    
    class outlier:
        abs_mag= []
        abs_mag_err = []
        period = []
        period_err = []
        log_period = []
        log_period_err = []
        obj = []
        diff_inerr = []


    """for i in range(len(objs)):
        if objs[i] == 'v396_cyg':
            outlier.abs_mag = np.append(outlier.abs_mag,abs_mags[i])
            outlier.abs_mag_err = np.append(outlier.abs_mag_err,abs_mags_errs[i])
            outlier.period = np.append(outlier.period,periods[i])
            outlier.period_err = np.append(outlier.period_err,periods_errs[i])
            outlier.obj = np.append(outlier.obj,objs[i])
            abs_mags = np.delete(abs_mags,i)
            periods = np.delete(periods,i)
            abs_mags_errs = np.delete(abs_mags_errs,i)
            periods_errs = np.delete(periods_errs,i)
            objs = np.delete(objs,i)
            break
        
    for i in range(len(objs)):
        if objs[i] == 'v1467_cyg':
            outlier.abs_mag = np.append(outlier.abs_mag,abs_mags[i])
            outlier.abs_mag_err = np.append(outlier.abs_mag_err,abs_mags_errs[i])
            outlier.period = np.append(outlier.period,periods[i])
            outlier.period_err = np.append(outlier.period_err,periods_errs[i])
            outlier.obj = np.append(outlier.obj,objs[i])
            abs_mags = np.delete(abs_mags,i)
            periods = np.delete(periods,i)
            abs_mags_errs = np.delete(abs_mags_errs,i)
            periods_errs = np.delete(periods_errs,i)
            objs = np.delete(objs,i)
            break
        
    for i in range(len(objs)):
        if objs[i] == 'v621_cyg':
            outlier.abs_mag = np.append(outlier.abs_mag,abs_mags[i])
            outlier.abs_mag_err = np.append(outlier.abs_mag_err,abs_mags_errs[i])
            outlier.period = np.append(outlier.period,periods[i])
            outlier.period_err = np.append(outlier.period_err,periods_errs[i])
            outlier.obj = np.append(outlier.obj,objs[i])
            abs_mags = np.delete(abs_mags,i)
            periods = np.delete(periods,i)
            abs_mags_errs = np.delete(abs_mags_errs,i)
            periods_errs = np.delete(periods_errs,i)
            objs = np.delete(objs,i)
            break
        
    for i in range(len(objs)):
        if objs[i] == 'v532_cyg':
            outlier.abs_mag = np.append(outlier.abs_mag,abs_mags[i])
            outlier.abs_mag_err = np.append(outlier.abs_mag_err,abs_mags_errs[i])
            outlier.period = np.append(outlier.period,periods[i])
            outlier.period_err = np.append(outlier.period_err,periods_errs[i])
            outlier.obj = np.append(outlier.obj,objs[i])
            abs_mags = np.delete(abs_mags,i)
            periods = np.delete(periods,i)
            abs_mags_errs = np.delete(abs_mags_errs,i)
            periods_errs = np.delete(periods_errs,i)
            objs = np.delete(objs,i)
            break
        
    for i in range(len(objs)):
        if objs[i] == 'kx_cyg':
            outlier.abs_mag = np.append(outlier.abs_mag,abs_mags[i])
            outlier.abs_mag_err = np.append(outlier.abs_mag_err,abs_mags_errs[i])
            outlier.period = np.append(outlier.period,periods[i])
            outlier.period_err = np.append(outlier.period_err,periods_errs[i])
            outlier.obj = np.append(outlier.obj,objs[i])
            abs_mags = np.delete(abs_mags,i)
            periods = np.delete(periods,i)
            abs_mags_errs = np.delete(abs_mags_errs,i)
            periods_errs = np.delete(periods_errs,i)
            objs = np.delete(objs,i)
            break
        
    for i in range(len(objs)):
        if objs[i] == 'sw_cas':
            outlier.abs_mag = np.append(outlier.abs_mag,abs_mags[i])
            outlier.abs_mag_err = np.append(outlier.abs_mag_err,abs_mags_errs[i])
            outlier.period = np.append(outlier.period,periods[i])
            outlier.period_err = np.append(outlier.period_err,periods_errs[i])
            outlier.obj = np.append(outlier.obj,objs[i])
            abs_mags = np.delete(abs_mags,i)
            periods = np.delete(periods,i)
            abs_mags_errs = np.delete(abs_mags_errs,i)
            periods_errs = np.delete(periods_errs,i)
            objs = np.delete(objs,i)
            break
        
    for i in range(len(objs)):
        if objs[i] == 'v609_cyg':
            outlier.abs_mag = np.append(outlier.abs_mag,abs_mags[i])
            outlier.abs_mag_err = np.append(outlier.abs_mag_err,abs_mags_errs[i])
            outlier.period = np.append(outlier.period,periods[i])
            outlier.period_err = np.append(outlier.period_err,periods_errs[i])
            outlier.obj = np.append(outlier.obj,objs[i])
            abs_mags = np.delete(abs_mags,i)
            periods = np.delete(periods,i)
            abs_mags_errs = np.delete(abs_mags_errs,i)
            periods_errs = np.delete(periods_errs,i)
            objs = np.delete(objs,i)
            break

    for i in range(len(objs)):
        if objs[i] == 'v538_cyg':
            outlier.abs_mag = np.append(outlier.abs_mag,abs_mags[i])
            outlier.abs_mag_err = np.append(outlier.abs_mag_err,abs_mags_errs[i])
            outlier.period = np.append(outlier.period,periods[i])
            outlier.period_err = np.append(outlier.period_err,periods_errs[i])
            outlier.obj = np.append(outlier.obj,objs[i])
            abs_mags = np.delete(abs_mags,i)
            periods = np.delete(periods,i)
            abs_mags_errs = np.delete(abs_mags_errs,i)
            periods_errs = np.delete(periods_errs,i)
            objs = np.delete(objs,i)
            break"""

    #outlier.log_period = np.log10(outlier.period)
    #outlier.log_period_err = outlier.period_err/(np.log(10)*outlier.period)
        
    log_periods = np.log10(periods)
    log_periods_errs = periods_errs/(np.log(10)*periods)
    
    initial_values = np.array([-3,1])
    
    print('initial Chi Squared = ', end='') # end='' to not start a new line - value printed by mean_absolute_deviation function
    initial_chi_squared = chi_squared(initial_values, fitting_model, periods, abs_mags, abs_mags_errs)
    
    fit = sp.optimize.minimize(chi_squared, # the function to minimize
                              initial_values, # where in 'parameter space' to start from
                              args=(fitting_model, periods, abs_mags, abs_mags_errs)) # model function and data to use
                             
    # Termination output message is fit.message - did the minimisation complete successfully?
    print(fit.message)
    
    print('minimised Chi Squared = {}'.format(fit.fun))
    chi_squared_min = chi_squared([fit.x[0], fit.x[1]], fitting_model, periods, abs_mags, abs_mags_errs)
    
    degrees_of_freedom = periods.size - fit.x.size # Make sure you understand why!
    print('DoF = {}'.format(degrees_of_freedom))
    
    reduced_chi_squared = chi_squared_min/degrees_of_freedom
    print('reduced chi^2 = {}'.format(reduced_chi_squared))
    
    P_value = sp.stats.chi2.sf(chi_squared_min, degrees_of_freedom)
    print('P(chi^2_min, DoF) = {}'.format(P_value))
    
    print('Slope: ' + str(fit.x[0]))
    
    smooth_y = fit.x[0] * log_periods + fit.x[1]
    
    diff = abs_mags - (fit.x[0]*log_periods + fit.x[1])
    
    diff_inerrs = diff/(abs_mags_errs)
    
    #outlier.diff_inerr = (outlier.abs_mag - (fit.x[0]*outlier.log_period + fit.x[1]))/outlier.abs_mag_err
    
    return objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, smooth_y, diff_inerrs
    
    


def jackknifing(abs_mags, abs_mags_errs, periods):
    
    jackknifed_slopes = []
    jackknifed_offsets = []
    
    for i in range(len(abs_mags)):
        jackknifed_abs_mags = np.delete(abs_mags,i)
        jackknifed_abs_mags_errs = np.delete(abs_mags_errs,i)
        jackknifed_periods = np.delete(periods,i)
        fit = sp.optimise.curvefit(chi_squared, np.array([-3,1]), args=(fitting_model, jackknifed_periods, jackknifed_abs_mags, jackknifed_abs_mags_errs))
        
        jackknifed_slopes = np.append(jackknifed_slopes,fit.x[0])
        jackknifed_offsets = np.append(jackknifed_offsets,fit.x[1])
        
    mean_slope=np.percentile(jackknifed_slopes,50)
    mean_offset=np.percentile(jackknifed_offsets,50)
    
    err_slope=np.std(jackknifed_periods)/np.sqrt(len(jackknifed_periods))
    err_offset=np.std(jackknifed_offsets)/np.sqrt(len(jackknifed_offsets))
    
    print("Slope: " + str(mean_slope) + " +/- " + str(err_slope))
    print("Offset: " + str(mean_offset) + " +/- " + str(err_offset))
    
    return np.array([mean_slope,err_slope]), np.array([mean_offset,err_offset])


def plot_pl(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, smooth_y, diff_inerrs):
    
    fig, axs = plt.subplots(2,1,height_ratios=(3,2))
    
    axs[1].set_xlabel(r"$log_{10}$(Period (days))")
    axs[0].set_ylabel("Absolute Magnitude")
    axs[0].set_xticklabels([])
    #axs[1].set_ylim([-5,5])
    #axs[1].set_yticks([-3,0,3],["-3","0","3"])
    axs[0].errorbar(log_periods, abs_mags, xerr=log_periods_errs, yerr=abs_mags_errs, linestyle=" ", marker='x', capsize=5, color="black")
    #axs[0].errorbar(outlier.log_period, outlier.abs_mag, xerr=outlier.log_period_err, yerr=outlier.abs_mag_err, linestyle=" ", marker='s', capsize=5, color="purple")
    axs[0].yaxis.set_inverted(True)
    axs[1].yaxis.set_inverted(True)
    axs[0].sharex(axs[1])
    axs[1].set_ylabel("Standard Errors")
    
    for i in range(len(objs)):
        axs[0].text(log_periods[i],abs_mags[i],objs[i],size=10)
        
    """for i in range(len(outlier.obj)):
        axs[0].text(outlier.log_period[i],outlier.abs_mag[i],outlier.obj[i],size=10,color="purple")"""
    
    axs[0].plot(log_periods, smooth_y, color='red', label='Our fit')
    
    gaia_y = -2.2 * log_periods - 2.05
    
    axs[0].plot(log_periods, gaia_y, color='green', label='Groenewegen (2018)')
    
    """fritz_y = -2.43 * smooth_x - 4.05
    
    plt.plot(smooth_x, fritz_y, color='blue', label='Fritz (2007)')"""
    
    axs[1].errorbar(log_periods,diff_inerrs,xerr=log_periods_errs,yerr=1,marker='x',linestyle=" ", capsize=5, color="black")
    
    #axs[1].errorbar(outlier.log_period,outlier.diff_inerr,xerr=outlier.log_period_err,yerr=1,marker='s',linestyle=" ", capsize=5, color="purple")
    
    axs[1].plot(log_periods,0*log_periods,color='black')
    #axs[1].plot(log_periods,(0*log_periods+3),color='gray')
    #axs[1].plot(log_periods,(0*log_periods-3),color='gray')
    
    axs[0].legend(loc='best')
    
    fig.subplots_adjust(hspace=0)
    
    plt.show()
    
    

def gaussian(x, xbar, sigma):
    return (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-xbar)**2/(2*sigma**2)))

def erf(x1, xbar, sigma):
    return sp.integrate.quad(lambda x: gaussian(x, xbar, sigma), -np.inf, x1)[0]
    
def chauvinet_criterion(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection):
    cycle = True
    
    outliers = []
    
    while cycle == True:
        
        objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, smooth_y, diff_inerrs = plot_calc(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection)
        
        maxindex = np.argmax(diff_inerrs)
        
        print(objs[maxindex])
        
        stddev = abs_mags_errs[maxindex]
        
        p_out = 0 - (erf((smooth_y[maxindex]+abs_mags[maxindex]),smooth_y[maxindex],stddev) - erf((smooth_y[maxindex]-abs_mags[maxindex]),smooth_y[maxindex],stddev))
        
        print(erf((smooth_y[maxindex]+abs_mags[maxindex]),smooth_y[maxindex],stddev))
        print(erf((smooth_y[maxindex]-abs_mags[maxindex]),smooth_y[maxindex],stddev))
        
        n_out = p_out * len(objs)
        print(n_out)
        
        if n_out < 0.5 :
            outliers = np.append(outliers,objs[maxindex])
            mags = np.delete(mags,maxindex)
            mags_errs = np.delete(mags_errs,maxindex)
            periods = np.delete(periods,maxindex)
            periods_errs = np.delete(periods_errs,maxindex)
            dist_pc = np.delete(dist_pc,maxindex)
            dist_pc_errs = np.delete(dist_pc_errs,maxindex)
            bvcorrection = np.delete(bvcorrection,maxindex)
            
        else:
            cycle = False
            
    return objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, smooth_y, diff_inerrs, outliers
        
    
mags, mags_errs, periods, periods_errs = get_periods()

#print(mags)
#print(mags_errs)

#dist_pc, dist_pc_errs = get_parallaxes()

dist_pc, dist_pc_errs, bvcorrection = read_parallaxes()

#objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, smooth_y, diff_inerrs, outliers = chauvinet_criterion(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection)

objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, smooth_y, diff_inerrs = plot_calc(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection)

plot_pl(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, smooth_y, diff_inerrs)



#print(erf(-1,-2,1))





