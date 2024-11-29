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
            break
        
        
    for i in range(len(objs)):
        if objs[i] == 'v407_cas':
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
    
    diff = abs_mags - (fit.x[0]*log_periods + fit.x[1])
    
    diff_inerrs = diff/(abs_mags_errs)
    
    #outlier.diff_inerr = (outlier.abs_mag - (fit.x[0]*outlier.log_period + fit.x[1]))/outlier.abs_mag_err
    
    return objs, periods, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier
    
def plot_calc_rm(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection):
    
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
        
    
    kill_list = np.array(["ak_cep","cd_cas","cg_cas","cy_cas","ry_cas","sz_cas","v1467_cyg","v342_cas","v538_cyg","v621_cyg","v609_cyg","vx_cyg"])
    
    for j in range(len(kill_list)):
        for i in range(len(objs)):
            if objs[i] == kill_list[j]:
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
    
    diff = abs_mags - (fit.x[0]*log_periods + fit.x[1])
    
    diff_inerrs = diff/(abs_mags_errs)
    
    #outlier.diff_inerr = (outlier.abs_mag - (fit.x[0]*outlier.log_period + fit.x[1]))/outlier.abs_mag_err
    
    return objs, periods, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier


def jackknifing(abs_mags, abs_mags_errs, periods):
    
    jackknifed_slopes = []
    jackknifed_offsets = []
    
    for i in range(len(abs_mags)):
        jackknifed_abs_mags = np.delete(abs_mags,i)
        jackknifed_abs_mags_errs = np.delete(abs_mags_errs,i)
        jackknifed_periods = np.delete(periods,i)
        initial_values = np.array([-3,1])
        fit = sp.optimize.minimize(chi_squared,
                                    initial_values,
                                    args=(fitting_model, jackknifed_periods, jackknifed_abs_mags, jackknifed_abs_mags_errs))
        
        jackknifed_slopes = np.append(jackknifed_slopes,fit.x[0])
        jackknifed_offsets = np.append(jackknifed_offsets,fit.x[1])
        
    mean_slope=np.mean(jackknifed_slopes)
    mean_offset=np.mean(jackknifed_offsets)
    
    err_slope=np.std(jackknifed_slopes)
    err_offset=np.std(jackknifed_offsets)
    
    print("Slope: " + str(mean_slope) + " +/- " + str(err_slope))
    print("Offset: " + str(mean_offset) + " +/- " + str(err_offset))
    
    return np.array([mean_slope,err_slope]), np.array([mean_offset,err_offset])



def jackgen(slope, offset, log_periods):
    
    jacky = slope[0] * log_periods + offset[0]
    
    return jacky



def plot_pl(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier, jacky):
    
    fig, axs = plt.subplots(2,1,height_ratios=(4,3))
    
    fig.set_size_inches(12,9)
    
    axs[1].set_xlabel(r"$log_{10}$(Period (days))")
    axs[0].set_ylabel("Absolute Magnitude")
    axs[0].set_xticklabels([])
    #axs[1].set_ylim([-5,5])
    #axs[1].set_yticks([-3,0,3],["-3","0","3"])
    axs[0].errorbar(log_periods, abs_mags, xerr=log_periods_errs, yerr=abs_mags_errs, linestyle=" ", marker='x', capsize=5, color="black")
    #axs[0].errorbar(outlier.log_period, outlier.abs_mag, xerr=outlier.log_period_err, yerr=outlier.abs_mag_err, linestyle=" ", marker='s', capsize=5, color="purple")       #Outlier
    axs[0].yaxis.set_inverted(True)
    axs[1].yaxis.set_inverted(True)
    axs[0].sharex(axs[1])
    axs[1].set_ylabel("Standard Errors")
    
    for i in range(len(objs)):
        axs[0].text(log_periods[i],abs_mags[i],objs[i],size=10)
        
    """for i in range(len(outlier.obj)):
        axs[0].text(outlier.log_period[i],outlier.abs_mag[i],outlier.obj[i],size=10,color="purple")  """               #Outlier
    
    axs[0].plot(log_periods, jacky, color='blue', label='Our fit')
    
    gaia_y = -2.2 * log_periods - 2.05
    
    axs[0].plot(log_periods, gaia_y, color='green', label='Groenewegen (2018)')
    
    """fritz_y = -2.43 * smooth_x - 4.05
    
    plt.plot(smooth_x, fritz_y, color='blue', label='Fritz (2007)')"""
    
    axs[1].errorbar(log_periods,diff_inerrs,xerr=log_periods_errs,yerr=1,marker='x',linestyle=" ", capsize=5, color="black")
    
    #axs[1].errorbar(outlier.log_period,outlier.diff_inerr,xerr=outlier.log_period_err,yerr=1,marker='s',linestyle=" ", capsize=5, color="purple")           #Outlier
    
    axs[1].plot(log_periods,0*log_periods,color='black')
    axs[1].plot(log_periods,(0*log_periods+1),color='gray')
    axs[1].plot(log_periods,(0*log_periods-1),color='gray')
    
    axs[0].legend(loc='best')
    
    fig.subplots_adjust(hspace=0)
    
    plt.savefig("plrelationship.png")
    plt.show()
    
def plot_pl_rm(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier, jacky):
    
    fig, axs = plt.subplots(2,1,height_ratios=(4,3))
    
    fig.set_size_inches(12,9)
    
    axs[1].set_xlabel(r"$log_{10}$(Period (days))")
    axs[0].set_ylabel("Absolute Magnitude")
    axs[0].set_xticklabels([])
    #axs[1].set_ylim([-5,5])
    #axs[1].set_yticks([-3,0,3],["-3","0","3"])
    axs[0].errorbar(log_periods, abs_mags, xerr=log_periods_errs, yerr=abs_mags_errs, linestyle=" ", marker='x', capsize=5, color="black")
    #axs[0].errorbar(outlier.log_period, outlier.abs_mag, xerr=outlier.log_period_err, yerr=outlier.abs_mag_err, linestyle=" ", marker='s', capsize=5, color="purple")       #Outlier
    axs[0].yaxis.set_inverted(True)
    axs[1].yaxis.set_inverted(True)
    axs[0].sharex(axs[1])
    axs[1].set_ylabel("Standard Errors")
    
    for i in range(len(objs)):
        axs[0].text(log_periods[i],abs_mags[i],objs[i],size=10)
        
    """for i in range(len(outlier.obj)):
        axs[0].text(outlier.log_period[i],outlier.abs_mag[i],outlier.obj[i],size=10,color="purple")  """               #Outlier
    
    axs[0].plot(log_periods, jacky, color='blue', label='Our fit')
    
    gaia_y = -2.2 * log_periods - 2.05
    
    axs[0].plot(log_periods, gaia_y, color='green', label='Groenewegen (2018)')
    
    """fritz_y = -2.43 * smooth_x - 4.05
    
    plt.plot(smooth_x, fritz_y, color='blue', label='Fritz (2007)')"""
    
    axs[1].errorbar(log_periods,diff_inerrs,xerr=log_periods_errs,yerr=1,marker='x',linestyle=" ", capsize=5, color="black")
    
    #axs[1].errorbar(outlier.log_period,outlier.diff_inerr,xerr=outlier.log_period_err,yerr=1,marker='s',linestyle=" ", capsize=5, color="purple")           #Outlier
    
    axs[1].plot(log_periods,0*log_periods,color='black')
    axs[1].plot(log_periods,(0*log_periods+1),color='gray')
    axs[1].plot(log_periods,(0*log_periods-1),color='gray')
    
    axs[0].legend(loc='best')
    
    fig.subplots_adjust(hspace=0)
    
    plt.savefig("plrelationship_sans_outliers.png")
    plt.show()

def aliasing(log_periods, abs_mags, abs_mags_errs, slope, offset):
    
    linslopes = np.linspace(-4,1,1001)
    linoffsets = np.linspace(-5,0,1001)
    
    slopechis = []
    offsetchis = []
    
    for linslope in linslopes:
        
        def linfunc(x, *param):
            return linslope*x + param[0]
        
        fit = sp.optimize.curve_fit(linfunc, log_periods, abs_mags, sigma=abs_mags_errs, absolute_sigma=True, p0=np.array([-3,1]), check_finite=True, maxfev=10**6)
        
        print(fit[0][0])
        
        minchi = chi_squared(fit[0][0],linfunc,log_periods,abs_mags,abs_mags_errs)
        
        slopechis = np.append(slopechis,minchi)
    
    
    for linoffset in linoffsets:
        
        def linfunc(x, *param):
            return param[0]*x + linoffset
        
        fit = sp.optimize.curve_fit(linfunc, log_periods, abs_mags, sigma=abs_mags_errs, absolute_sigma=True, p0=np.array([-3,1]), check_finite=True, maxfev=10**6)
        
        minchi = chi_squared(fit[0][0],linfunc,log_periods,abs_mags,abs_mags_errs)
        
        offsetchis = np.append(offsetchis,minchi)
        
    plt.plot(linslopes,slopechis)
    plt.xlabel("Slope")
    plt.ylabel("Minimised Chi-Squared")
    plt.show()
    plt.plot(linoffsets,offsetchis)
    plt.xlabel("Offset")
    plt.ylabel("Minimised Chi-Squared")
    plt.show()
    
    def linfunc(x, *params):
            return params[0][0]*x + params[0][1]
        
    chi_squared_min = chi_squared([slope[0],offset[0]], linfunc, log_periods, abs_mags, abs_mags_errs)
    
    extent = 3.5        #Standard Errors
    n_points = 150      #Mesh Density
    
    p0_range = extent * slope[1]
    p1_range = extent * offset[1]
    
    #Generate grid and data
    p0_axis = np.linspace(slope[0]-p0_range, slope[0]+p0_range, num=n_points)
    p1_axis = np.linspace(offset[0]-p1_range, offset[0]+p1_range, num=n_points)
    plot_data = np.zeros((n_points,n_points))
    
    for j, p1_val in enumerate(p1_axis): 
        for i, p0_val in enumerate(p0_axis): # Nested loops for 'clarity'...
            plot_data[j][i] = chi_squared([p0_val, p1_val], # function evaluated n_points*n_points times!
                                        linfunc, 
                                        log_periods, 
                                        abs_mags, 
                                        abs_mags_errs)
            
    plt.figure(figsize=(10,8))
    im = plt.imshow(plot_data, # grid of chi-squared values
                    extent=(p0_axis[0], p0_axis[-1], # 'x' range
                            p1_axis[0], p1_axis[-1]), # 'y' range
                    origin='lower', aspect='auto')

    plt.xlim(slope[0]-p0_range, slope[0]+p0_range) # axis ranges
    plt.ylim(offset[0]-p1_range, offset[0]+p1_range)

    plt.xlabel('Slope') # Axis labels
    plt.ylabel('Offset')

    cbar=plt.colorbar(im, orientation='vertical') # Colorbar and label
    cbar.set_label('$\chi^2$', fontsize=12)

    # Add in best fit point and dashed lines to the axes
    plt.plot(slope[0], offset[0], 'wo') 
    plt.plot((slope[0], slope[0]), (p1_axis[0], offset[0]), # vertical line
            linestyle='--', color='w')
    plt.plot((p0_axis[0], slope[0]), (offset[0], offset[0]), # horizontal line
            linestyle='--', color='w')
    plt.savefig("pl_chisquared_heatmap.png")
    plt.show()
    
    X, Y = np.meshgrid(p0_axis, p1_axis, indexing='xy')
    contour_data = plot_data - chi_squared_min

    levels = [1, 4, 9] # Contour levels in delta chi-squared of 1, 4 & 9 correspond to 1, 2 & 3 standard errors
    plt.figure(figsize=(12,8))

    #Plot and label contours: (comment out labelling to remove text over contours)
    contour_plot = plt.contour(X, Y, contour_data, levels=levels, colors='b', origin='lower')
    plt.clabel(contour_plot, levels, fontsize=12, inline=1, fmt=r'$\chi^2 = \chi^2_{min}+%1.0f$') 

    plt.xlabel('Slope') 
    plt.ylabel('Offset')

    import matplotlib.ticker as ticker # Allows you to modify the tick markers to 
    xtick_spacing = p0_range/4         # assess the errors from the chi-squared 
    ytick_spacing = p1_range/4         # contour plots - set as appropriate
    
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
    ax.xaxis.set_major_formatter('{x:.2f}') # 2 decimal places - may or may not be appropriate!

    ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
    ax.yaxis.set_major_formatter('{x:.2f}') # 2 decimal places - may or may not be appropriate!

    # Add in best fit point and dashed lines to the axes
    plt.plot(slope[0], offset[0], 'ro') 
    plt.plot((slope[0], slope[0]), (p1_axis[0], offset[0]), linestyle='--', color='r')
    plt.plot((p0_axis[0], slope[0]), (offset[0], offset[0]), linestyle='--', color='r')
    plt.savefig("pl_chisquared_contour.png")
    plt.show()
    
def aliasing_rm(log_periods, abs_mags, abs_mags_errs, slope, offset):
    
    linslopes = np.linspace(-4,1,1001)
    linoffsets = np.linspace(-5,0,1001)
    
    slopechis = []
    offsetchis = []
    
    for linslope in linslopes:
        
        def linfunc(x, *param):
            return linslope*x + param[0]
        
        fit = sp.optimize.curve_fit(linfunc, log_periods, abs_mags, sigma=abs_mags_errs, absolute_sigma=True, p0=np.array([-3,1]), check_finite=True, maxfev=10**6)
        
        print(fit[0][0])
        
        minchi = chi_squared(fit[0][0],linfunc,log_periods,abs_mags,abs_mags_errs)
        
        slopechis = np.append(slopechis,minchi)
    
    
    for linoffset in linoffsets:
        
        def linfunc(x, *param):
            return param[0]*x + linoffset
        
        fit = sp.optimize.curve_fit(linfunc, log_periods, abs_mags, sigma=abs_mags_errs, absolute_sigma=True, p0=np.array([-3,1]), check_finite=True, maxfev=10**6)
        
        minchi = chi_squared(fit[0][0],linfunc,log_periods,abs_mags,abs_mags_errs)
        
        offsetchis = np.append(offsetchis,minchi)
        
    plt.plot(linslopes,slopechis)
    plt.xlabel("Slope")
    plt.ylabel("Minimised Chi-Squared")
    plt.show()
    plt.plot(linoffsets,offsetchis)
    plt.xlabel("Offset")
    plt.ylabel("Minimised Chi-Squared")
    plt.show()
    
    def linfunc(x, *params):
            return params[0][0]*x + params[0][1]
        
    chi_squared_min = chi_squared([slope[0],offset[0]], linfunc, log_periods, abs_mags, abs_mags_errs)
    
    extent = 3.5        #Standard Errors
    n_points = 150      #Mesh Density
    
    p0_range = extent * slope[1]
    p1_range = extent * offset[1]
    
    #Generate grid and data
    p0_axis = np.linspace(slope[0]-p0_range, slope[0]+p0_range, num=n_points)
    p1_axis = np.linspace(offset[0]-p1_range, offset[0]+p1_range, num=n_points)
    plot_data = np.zeros((n_points,n_points))
    
    for j, p1_val in enumerate(p1_axis): 
        for i, p0_val in enumerate(p0_axis): # Nested loops for 'clarity'...
            plot_data[j][i] = chi_squared([p0_val, p1_val], # function evaluated n_points*n_points times!
                                        linfunc, 
                                        log_periods, 
                                        abs_mags, 
                                        abs_mags_errs)
            
    plt.figure(figsize=(10,8))
    im = plt.imshow(plot_data, # grid of chi-squared values
                    extent=(p0_axis[0], p0_axis[-1], # 'x' range
                            p1_axis[0], p1_axis[-1]), # 'y' range
                    origin='lower', aspect='auto')

    plt.xlim(slope[0]-p0_range, slope[0]+p0_range) # axis ranges
    plt.ylim(offset[0]-p1_range, offset[0]+p1_range)

    plt.xlabel('Slope') # Axis labels
    plt.ylabel('Offset')

    cbar=plt.colorbar(im, orientation='vertical') # Colorbar and label
    cbar.set_label('$\chi^2$', fontsize=12)

    # Add in best fit point and dashed lines to the axes
    plt.plot(slope[0], offset[0], 'wo') 
    plt.plot((slope[0], slope[0]), (p1_axis[0], offset[0]), # vertical line
            linestyle='--', color='w')
    plt.plot((p0_axis[0], slope[0]), (offset[0], offset[0]), # horizontal line
            linestyle='--', color='w')
    plt.savefig("pl_chisquared_heatmap_sans_outliers.png")
    plt.show()
    
    X, Y = np.meshgrid(p0_axis, p1_axis, indexing='xy')
    contour_data = plot_data - chi_squared_min

    levels = [1, 4, 9] # Contour levels in delta chi-squared of 1, 4 & 9 correspond to 1, 2 & 3 standard errors
    plt.figure(figsize=(12,8))

    #Plot and label contours: (comment out labelling to remove text over contours)
    contour_plot = plt.contour(X, Y, contour_data, levels=levels, colors='b', origin='lower')
    plt.clabel(contour_plot, levels, fontsize=12, inline=1, fmt=r'$\chi^2 = \chi^2_{min}+%1.0f$') 

    plt.xlabel('Slope') 
    plt.ylabel('Offset')

    import matplotlib.ticker as ticker # Allows you to modify the tick markers to 
    xtick_spacing = p0_range/4         # assess the errors from the chi-squared 
    ytick_spacing = p1_range/4         # contour plots - set as appropriate
    
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
    ax.xaxis.set_major_formatter('{x:.2f}') # 2 decimal places - may or may not be appropriate!

    ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
    ax.yaxis.set_major_formatter('{x:.2f}') # 2 decimal places - may or may not be appropriate!

    # Add in best fit point and dashed lines to the axes
    plt.plot(slope[0], offset[0], 'ro') 
    plt.plot((slope[0], slope[0]), (p1_axis[0], offset[0]), linestyle='--', color='r')
    plt.plot((p0_axis[0], slope[0]), (offset[0], offset[0]), linestyle='--', color='r')
    plt.savefig("pl_chisquared_contour_sans_outliers.png")
    plt.show()    


def pl_gen():    
    mags, mags_errs, periods, periods_errs = get_periods()
    #dist_pc, dist_pc_errs = get_parallaxes()
    dist_pc, dist_pc_errs, bvcorrection = read_parallaxes()
    objs, periods, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier = plot_calc(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection)
    slope_array, offset_array = jackknifing(abs_mags, abs_mags_errs, periods)
    jacky = jackgen(slope_array, offset_array, log_periods)
    plot_pl(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier, jacky)
    aliasing(log_periods, abs_mags, abs_mags_errs, slope_array, offset_array)

def pl_gen_sans_outliers():
    mags, mags_errs, periods, periods_errs = get_periods()

    #dist_pc, dist_pc_errs = get_parallaxes()

    dist_pc, dist_pc_errs, bvcorrection = read_parallaxes()

    objs, periods, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier = plot_calc_rm(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection)

    slope_array, offset_array = jackknifing(abs_mags, abs_mags_errs, periods)

    jacky = jackgen(slope_array, offset_array, log_periods)

    plot_pl_rm(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier, jacky)

    aliasing_rm(log_periods, abs_mags, abs_mags_errs, slope_array, offset_array)


pl_gen()
pl_gen_sans_outliers()