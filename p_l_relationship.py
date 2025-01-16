from imports import *
from data_handling import get_bv
plt.rcParams.update({'font.size': 22})
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)


def get_periods():
    
    df = pd.read_csv("periodmagdata.txt", delimiter="x", header=None)
    
    mags = df[0].to_numpy()
    mags_errs = df[1].to_numpy()
    periods = df[2].to_numpy()
    str_periods_errs = df[3].to_numpy()

    periods_errs=[[],[]]

    for i in str_periods_errs:
        arr=literal_eval(i)
        negative_error=abs(arr[0])
        positive_error=abs(arr[1])

        periods_errs[0].append(negative_error)
        periods_errs[1].append(positive_error)

    periods_errs=np.array(periods_errs)
    
    return mags, mags_errs, periods, periods_errs

def get_bvcorrection():
    
    objs=next(os.walk('.'))[1]
    objs = objs[3:-1]
    
    Vmag = []
    Bmag = []
    
    Verr = []
    Berr = []
    
    for obj in objs:
        newB, newV = get_bv(obj)
        
        newVerr = np.std(newV)
        newBerr = np.std(newB)
        
        newB = np.mean(newB)
        newV = np.mean(newV)
        
        Vmag = np.append(Vmag, newV)
        Bmag = np.append(Bmag, newB)
        
        Verr = np.append(Verr,newVerr)
        Berr = np.append(Berr,newBerr)
        
    bvcorrection = Bmag - Vmag
    
    df = pd.read_csv("intrinsic.txt",delimiter=" ")
    spectraltypedata = df[["Type","B-V"]].to_numpy()
    
    Simbad.add_votable_fields("sptype")
    
    for j in range(len(objs)):
        
        obj = objs[j]
        
        result = Simbad.query_object(obj)
        
        for i in range(len(spectraltypedata)):
            
            if spectraltypedata[i][0] == result['SP_TYPE'][0]:
                
                bvcorrection[j] = bvcorrection[j] - spectraltypedata[i][1]
    
    bverror = np.sqrt(Verr**2 + Berr**2)

    print(f'BV Correction: {bvcorrection}')
    
    return bvcorrection, bverror

def get_parallaxes():
    
    objs=next(os.walk('.'))[1]
    objs = objs[3:-1]
    
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
    
    bvcorrection, bverror = get_bvcorrection()
    
    return dist_pc, dist_pc_errs, bvcorrection, bverror



def fitting_model(x, params):
    return params[0]*np.log10(x) + params[1]

def chi_squared(model_params, model, x_data, y_data, y_err):
    chi_sq = np.sum(((y_data - model(x_data, model_params))/y_err)**2)
    #print(chi_sq)
    return chi_sq

def mean_absolute_deviation(model_params, model, x_data, y_data, npoints):
    mad = np.sum(np.abs(y_data - model(x_data, model_params)))/npoints
    #print(mad)
    return mad



def plot_calc(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection, bverror):
    
    objs=next(os.walk('.'))[1]
    objs = np.array(objs[3:-1])
    
    #print(bvcorrection)
    
    abs_mags = (mags - 3.14 * bvcorrection) - 5*np.log10(dist_pc) + 5
    print(mags_errs)
    print(bverror)
    print(5*(dist_pc_errs/(np.log(10)*dist_pc)))
    abs_mags_errs = np.sqrt(((mags_errs**2 + (3.14 * bverror)**2)) + (5*(dist_pc_errs/(np.log(10)*dist_pc)))**2)
    
    sort_index = np.argsort(periods)
    
    abs_mags = abs_mags[sort_index]
    periods = periods[sort_index]
    abs_mags_errs = abs_mags_errs[sort_index]
    periods_errs = [periods_errs[0][sort_index],periods_errs[1][sort_index]]
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
        
    log_periods = np.log10(periods)
    log_periods_errs = periods_errs/(np.log(10)*periods)
    
    initial_values = np.array([-2.43,-4.05])

    print(f'Correlation Coefficient: {np.corrcoef(log_periods, abs_mags)[0, 1]}')
    print(f'Pearson Correlation Coefficient: {sp.stats.pearsonr(log_periods, abs_mags)}')
    
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
    
def plot_calc_rm(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection, bverror):
    
    objs=next(os.walk('.'))[1]
    objs = np.array(objs[3:-1])
    
    #print(bvcorrection)
    
    abs_mags = (mags - 3.1* bvcorrection) - 5*np.log10(dist_pc) + 5
    
    abs_mags_errs = np.sqrt(((mags_errs**2 + (3.1 * bverror)**2)) + (5*(dist_pc_errs/(np.log(10)*dist_pc)))**2)
    
    sort_index = np.argsort(periods)
    
    abs_mags = abs_mags[sort_index]
    periods = periods[sort_index]
    abs_mags_errs = abs_mags_errs[sort_index]
    periods_errs = np.array([periods_errs[0][sort_index],periods_errs[1][sort_index]])
    objs = objs[sort_index]

    print(np.shape(periods_errs))

    class outlier:
        abs_mag= []
        abs_mag_err = []
        period = []
        period_err = []
        log_period = []
        log_period_err = []
        obj = []
        diff_inerr = []
        
    
    kill_list = np.array(['v1467_cyg',"v396_cyg"])
    
    for j in range(len(kill_list)):
        for i in range(len(objs)):
            if objs[i] == kill_list[j]:

                print(np.shape(periods_errs))

                outlier.abs_mag = np.append(outlier.abs_mag,abs_mags[i])
                outlier.abs_mag_err = np.append(outlier.abs_mag_err,abs_mags_errs[i])
                outlier.period = np.append(outlier.period,periods[i])
                outlier.period_err = np.append(outlier.period_err,[periods_errs[0][i],periods_errs[1][i]])
                outlier.obj = np.append(outlier.obj,objs[i])
                abs_mags = np.delete(abs_mags,i)
                periods = np.delete(periods,i)
                abs_mags_errs = np.delete(abs_mags_errs,i)
                top_errs = np.delete(periods_errs[0],i)
                bottom_errs = np.delete(periods_errs[1],i)
                periods_errs=[top_errs,bottom_errs]
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
        initial_values = np.array([-2.43,-4.05])
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
    chi_squared_min = chi_squared([fit.x[0], fit.x[1]], fitting_model, periods, abs_mags, abs_mags_errs)
    degrees_of_freedom = len(abs_mags) -2
    reduced_chi_squared = chi_squared_min/degrees_of_freedom
    
    sum_squared_errors = np.sum((abs_mags - (mean_slope*np.log10(periods) + mean_offset))**2)
    total_sum_squared = np.sum(abs_mags**2)
    
    r_squared = 1 - sum_squared_errors/total_sum_squared
    
    with open("PLOutputs.txt","w") as f:
        print("## All points fitting ##", file=f)
        print("Slope: " + str(mean_slope) + " +/- " + str(err_slope), file=f)
        print("Offset: " + str(mean_offset) + " +/- " + str(err_offset), file=f)
        print("Minimised Chi-Squared: " + str(chi_squared_min), file=f)
        print("Degrees of Freedom: " + str(degrees_of_freedom), file=f)
        print("Reduced Chi-Squared: " + str(reduced_chi_squared), file=f)
        print("R-squared for fit: " + str(r_squared), file=f)
        print(" ", file=f)
        f.close()
        
    return np.array([mean_slope,err_slope]), np.array([mean_offset,err_offset])

def jackknifing_a(abs_mags, abs_mags_errs, periods):
    
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
    chi_squared_min = chi_squared([fit.x[0], fit.x[1]], fitting_model, periods, abs_mags, abs_mags_errs)
    degrees_of_freedom = len(abs_mags) -2
    reduced_chi_squared = chi_squared_min/degrees_of_freedom
    
    sum_squared_errors = np.sum((abs_mags - (mean_slope*np.log10(periods) + mean_offset))**2)
    total_sum_squared = np.sum(abs_mags**2)
    
    r_squared = 1 - sum_squared_errors/total_sum_squared
    
    with open("PLOutputs.txt","a") as f:
        print("## Fitting parameters sans overfitting points ##", file=f)
        print("Slope: " + str(mean_slope) + " +/- " + str(err_slope), file=f)
        print("Offset: " + str(mean_offset) + " +/- " + str(err_offset), file=f)
        print("Minimised Chi-Squared: " + str(chi_squared_min), file=f)
        print("Degrees of Freedom: " + str(degrees_of_freedom), file=f)
        print("Reduced Chi-Squared: " + str(reduced_chi_squared), file=f)
        print("R-squared for fit: " + str(r_squared), file=f)
        print(" ", file=f)
        f.close()
        
    return np.array([mean_slope,err_slope]), np.array([mean_offset,err_offset])

def jackgen(slope, offset, log_periods):
    
    jacky = slope[0] * log_periods + offset[0]
    
    return jacky



def plot_pl(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier, jacky):
    
    fig, axs = plt.subplots(2,1,height_ratios=(3,1))
    
    fig.set_size_inches(12,9)
    
    axs[1].set_xlabel(r"$log_{10}$(Period (days))")
    axs[0].set_ylabel("Absolute Magnitude")
    axs[0].set_xticklabels([])
    #axs[1].set_ylim([-5,5])
    #axs[1].set_yticks([-3,0,3],["-3","0","3"])
    axs[0].errorbar(log_periods, abs_mags, xerr=log_periods_errs, yerr=abs_mags_errs, linestyle=" ", marker='x', capsize=5, color="black")
    #axs[0].errorbar(outlier.log_period, outlier.abs_mag, xerr=outlier.log_period_err, yerr=outlier.abs_mag_err, linestyle=" ", marker='s', capsize=5, color="lime")       #Outlier
    axs[0].yaxis.set_inverted(True)
    axs[1].yaxis.set_inverted(True)
    axs[0].sharex(axs[1])
    axs[1].set_ylabel("Norm. Residuals")

    slope, intercept, r, p, se = sp.stats.linregress(log_periods, abs_mags)

    axs[0].plot(log_periods,slope*log_periods+intercept,c='lime')

    axs[0].plot(log_periods, jacky, color='red', label='Our fit')
    
    """gaia_y = -2.2 * log_periods - 2.05
    
    axs[0].plot(log_periods, gaia_y, color='green', label='Groenewegen (2018)')"""
    
    fritz_y = -2.43 * log_periods - 4.05
    
    axs[0].plot(log_periods, fritz_y, color='blue', label='Fritz (2007)')
    
    axs[1].errorbar(log_periods,diff_inerrs,xerr=log_periods_errs,yerr=1,marker='x',linestyle=" ", capsize=5, color="r",ecolor='k')
    
    #axs[1].errorbar(outlier.log_period,outlier.diff_inerr,xerr=outlier.log_period_err,yerr=1,marker='s',linestyle=" ", capsize=5, color="purple")           #Outlier
    
    axs[1].axhline(0,color='black',linestyle='dashed')
    #axs[1].plot(log_periods,(0*log_periods+1),color='gray')
    #axs[1].plot(log_periods,(0*log_periods-1),color='gray')
    
    #axs[0].legend(loc='best')

    fig.subplots_adjust(right=0.8)
        
    reshist=fig.add_axes([0.8,0.11,0.1925,0.1925])

    reshist.hist(diff_inerrs,color='r',ec='r',orientation='horizontal')
    reshist.set_yticks([])
    
    fig.subplots_adjust(hspace=0)
    plt.savefig("plrelationship.png")
    plt.show()
    print(objs)
    print(len(objs))
    
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
    
    axs[0].plot(log_periods, jacky, color='red', label='Our fit')
    
    """gaia_y = -2.2 * log_periods - 2.05
    
    axs[0].plot(log_periods, gaia_y, color='green', label='Groenewegen (2018)')"""
    
    fritz_y = -2.43 * log_periods - 4.05
    
    axs[0].plot(log_periods, fritz_y, color='blue', label='Fritz (2007)')
    
    axs[1].errorbar(log_periods,diff_inerrs,xerr=log_periods_errs,yerr=1,marker='x',linestyle=" ", capsize=5, color="black")
    
    #axs[1].errorbar(outlier.log_period,outlier.diff_inerr,xerr=outlier.log_period_err,yerr=1,marker='s',linestyle=" ", capsize=5, color="purple")           #Outlier
    
    axs[1].plot(log_periods,0*log_periods,color='black')
    axs[1].plot(log_periods,(0*log_periods+1),color='gray')
    axs[1].plot(log_periods,(0*log_periods-1),color='gray')
    
    axs[0].legend(loc='best')
    
    fig.subplots_adjust(hspace=0)
    plt.savefig("plrelationship_sans_outliers.png")

def plot_pl_odr(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, beta):
    
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
    
    axs[0].plot(log_periods, beta[0]*log_periods+beta[1], color='r', label='Our fit')
    
    """gaia_y = -2.2 * log_periods - 2.05
    
    axs[0].plot(log_periods, gaia_y, color='green', label='Groenewegen (2018)')"""
    
    fritz_y = -2.43 * log_periods - 4.05
    
    axs[0].plot(log_periods, fritz_y, color='blue', label='Fritz (2007)')
    
    axs[1].errorbar(log_periods,diff_inerrs,xerr=log_periods_errs,yerr=1,marker='x',linestyle=" ", capsize=5, color="black")
    
    #axs[1].errorbar(outlier.log_period,outlier.diff_inerr,xerr=outlier.log_period_err,yerr=1,marker='s',linestyle=" ", capsize=5, color="purple")           #Outlier
    
    axs[1].plot(log_periods,0*log_periods,color='black')
    axs[1].plot(log_periods,(0*log_periods+1),color='gray')
    axs[1].plot(log_periods,(0*log_periods-1),color='gray')
    
    axs[0].legend(loc='best')
    
    fig.subplots_adjust(hspace=0)
    plt.savefig("plrelationship_odr.png")
    print(objs)
    print(len(objs))
    



def aliasing(log_periods, abs_mags, abs_mags_errs, slope, offset):
    
    linslopes = np.linspace(-4,1,1001)
    linoffsets = np.linspace(-5,0,1001)
    
    slopechis = []
    offsetchis = []
    
    for linslope in linslopes:
        
        def linfunc(x, *param):
            return linslope*x + param[0]
        
        fit = sp.optimize.curve_fit(linfunc, log_periods, abs_mags, sigma=abs_mags_errs, absolute_sigma=True, p0=np.array([-3,1]), check_finite=True, maxfev=10**6)
        
        #print(fit[0][0])
        
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
    plt.tight_layout()
    
    plt.plot(linoffsets,offsetchis)
    plt.xlabel("Offset")
    plt.ylabel("Minimised Chi-Squared")
    plt.tight_layout()
    
    
    def linfunc(x, *params):
            return params[0][0]*x + params[0][1]
        
    chi_squared_min = chi_squared([slope[0],offset[0]], linfunc, log_periods, abs_mags, abs_mags_errs)
    
    extent = 3        #Standard Errors
    n_points = 1000      #Mesh Density
    
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
    plt.tight_layout()
    
    
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
    plt.tight_layout()
    
    
def aliasing_rm(log_periods, abs_mags, abs_mags_errs, slope, offset):
    
    linslopes = np.linspace(-4,1,1001)
    linoffsets = np.linspace(-5,0,1001)
    
    slopechis = []
    offsetchis = []
    
    for linslope in linslopes:
        
        def linfunc(x, *param):
            return linslope*x + param[0]
        
        fit = sp.optimize.curve_fit(linfunc, log_periods, abs_mags, sigma=abs_mags_errs, absolute_sigma=True, p0=np.array([-3,1]), check_finite=True, maxfev=10**6)
        
        #print(fit[0][0])
        
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
    plt.tight_layout()
    
    plt.plot(linoffsets,offsetchis)
    plt.xlabel("Offset")
    plt.ylabel("Minimised Chi-Squared")
    plt.tight_layout()
    
    
    def linfunc(x, *params):
            return params[0][0]*x + params[0][1]
        
    chi_squared_min = chi_squared([slope[0],offset[0]], linfunc, log_periods, abs_mags, abs_mags_errs)
    
    extent = 3        #Standard Errors
    n_points = 1000      #Mesh Density
    
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
    plt.tight_layout()
    
    
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
    plt.tight_layout()
        

def aliasing_odr(log_periods, abs_mags, log_periods_errs, abs_mags_errs, slope, offset):
    def linfunc(p, x):
            m, c = p
            return m*x + c
    
    lin_model = sp.odr.Model(linfunc)
    
    data = sp.odr.RealData(log_periods, abs_mags, sx=log_periods_errs, sy=abs_mags_errs)
    
    odr = sp.odr.ODR(data, lin_model, beta0=[slope[0],offset[0]])
    
    out = odr.run()
    
    sum_square_error_min = out.res_var
    print("Min Res Var: " + str(sum_square_error_min))
    
    
    """testoffset = np.linspace(offset[0]-10,offset[0]+10,101)
    testminresvar = []
    for i in range(101):
        odr = sp.odr.ODR(data, lin_model, beta0=[slope[0],testoffset[i]],maxit=0)
        out = odr.run()
        testminresvar = np.append(testminresvar,out.res_var - sum_square_error_min)
        
    plt.plot(testoffset,testminresvar)
    """
    
    extent = 3        #Standard Errors
    n_points = 1000     #Mesh Density
    
    p0_range = extent * slope[1]
    p1_range = extent * offset[1]
    
    #Generate grid and data
    p0_axis = np.linspace(slope[0]-p0_range, slope[0]+p0_range, num=n_points)
    p1_axis = np.linspace(offset[0]-p1_range, offset[0]+p1_range, num=n_points)
    plot_data = np.zeros((n_points,n_points))
    
    for j, p1_val in enumerate(p1_axis): 
        for i, p0_val in enumerate(p0_axis): # Nested loops for 'clarity'...
            newodr = sp.odr.ODR(data, lin_model, beta0=[p0_val,p1_val], overwrite=True, maxit=0)
            out = newodr.run()
            plot_data[j][i] = out.res_var
            
    plt.figure(figsize=(10,8))
    im = plt.imshow(plot_data, # grid of chi-squared values
                    extent=(p0_axis[0], p0_axis[-1], # 'x' range
                            p1_axis[0], p1_axis[-1]), # 'y' range
                    origin='lower', aspect='auto',
                    cmap="YlOrBr_r")

    plt.xlim(slope[0]-p0_range, slope[0]+p0_range) # axis ranges
    plt.ylim(offset[0]-p1_range, offset[0]+p1_range)

    plt.xlabel('Slope') # Axis labels
    plt.ylabel('Offset')

    cbar=plt.colorbar(im, orientation='vertical') # Colorbar and label
    cbar.set_label('Residual Variance', fontsize=20)
    
    import matplotlib.ticker as ticker # Allows you to modify the tick markers to 
    xtick_spacing = p0_range/4         # assess the errors from the chi-squared 
    ytick_spacing = p1_range/4         # contour plots - set as appropriate
    
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
    ax.xaxis.set_major_formatter('{x:.2f}') # 2 decimal places - may or may not be appropriate!

    ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
    ax.yaxis.set_major_formatter('{x:.2f}') # 2 decimal places - may or may not be appropriate!

    # Add in best fit point and dashed lines to the axes
    plt.plot(slope[0], offset[0], color='White', marker='o') 
    plt.plot((slope[0], slope[0]), (p1_axis[0], offset[0]), # vertical line
            linestyle='--', color='White')
    plt.plot((p0_axis[0], slope[0]), (offset[0], offset[0]), # horizontal line
            linestyle='--', color='White')
    plt.savefig("pl_odr_sse_heatmap.png")
    plt.tight_layout()
    
    X, Y = np.meshgrid(p0_axis, p1_axis, indexing='xy')
    contour_data = plot_data - sum_square_error_min

    levels = [1, 4, 9] # Contour levels in delta chi-squared of 1, 4 & 9 correspond to 1, 2 & 3 standard errors
    plt.figure(figsize=(10,8))

    #Plot and label contours: (comment out labelling to remove text over contours)
    contour_plot = plt.contour(X, Y, contour_data, levels=levels, colors='Black', origin='lower')
    plt.clabel(contour_plot, levels, fontsize=18, inline=1, fmt=r'$RV = RV_{min}+%1.0f$') 

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
    plt.savefig("pl_odr_sse_contour.png")
    plt.tight_layout()



def orthogonal_distance_regression(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection, bverror):
    
    objs=next(os.walk('.'))[1]
    objs = np.array(objs[3:-1])
    
    #print(bvcorrection)
    
    abs_mags = (mags - 3.14 * bvcorrection) - 5*np.log10(dist_pc) + 5
    
    abs_mags_errs = np.sqrt(((mags_errs**2 + (3.14 * bverror)**2)) + (5*(dist_pc_errs/(np.log(10)*dist_pc)))**2)
    
    sort_index = np.argsort(periods)
    
    abs_mags = abs_mags[sort_index]
    periods = periods[sort_index]
    abs_mags_errs = abs_mags_errs[sort_index]
    sorted_periods_errs = np.array([periods_errs[0][sort_index],periods_errs[1][sort_index]])
    periods_errs=[]
    for i in range(len(sorted_periods_errs[0])):
        if sorted_periods_errs[0][i]>sorted_periods_errs[1][i]:
            periods_errs.append(sorted_periods_errs[0][i])
        else:
            periods_errs.append(sorted_periods_errs[1][i])

    objs = objs[sort_index]
    
    def log_func(p, x):
        m, c = p
        return m*np.log10(x) + c
    
    log_model = sp.odr.Model(log_func)
    
    data = sp.odr.RealData(periods, abs_mags, sx=periods_errs, sy=abs_mags_errs)
    
    odr = sp.odr.ODR(data, log_model, beta0=[-2.17175,3.57509])
    
    out = odr.run()
    
    out.pprint()
    
    log_periods = np.log10(periods)
    log_periods_errs = periods_errs/(np.log(10)*periods)

    diff = abs_mags - (out.beta[0]*log_periods + out.beta[1])
    
    diff_inerrs = diff/(abs_mags_errs)
    
    degrees_of_freedom = len(objs) - len(out.beta)
    
    sum_square_error = np.sum((abs_mags - (out.beta[0]*np.log10(periods) + out.beta[1]))**2)
    total_sum_square = np.sum(abs_mags**2)
    
    r_squared = 1 - sum_square_error/total_sum_square
    
    with open("PLOutputs.txt","a") as f:
        print("## Orthogonal Distance Regression Values ##", file=f)
        print("Slope: " + str(out.beta[0]) + " +/- " + str(out.sd_beta[0]), file=f)
        print("Offset: " + str(out.beta[1]) + " +/- " + str(out.sd_beta[1]), file=f)
        print("Minimised Sum of Squares of errors: " + str(out.sum_square), file=f)
        print("Degrees of Freedom: " + str(degrees_of_freedom), file=f)
        print("Residual Variance: " + str(out.sum_square/degrees_of_freedom), file=f)
        print("R-squared for fit: " + str(r_squared), file=f)
        f.close()
    
    return objs, periods, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, out.beta, np.array([out.beta[0],out.sd_beta[0]]), np.array([out.beta[1],out.sd_beta[1]])



def pl_gen():    
    mags, mags_errs, periods, periods_errs = get_periods()
    #dist_pc, dist_pc_errs = get_parallaxes()
    dist_pc, dist_pc_errs, bvcorrection, bverror = read_parallaxes()
    objs, periods, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier = plot_calc(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection, bverror)
    slope_array, offset_array = jackknifing(abs_mags, abs_mags_errs, periods)
    jacky = jackgen(slope_array, offset_array, log_periods)
    plot_pl(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier, jacky)
    aliasing(log_periods, abs_mags, abs_mags_errs, slope_array, offset_array)

def pl_gen_sans_outliers():
    mags, mags_errs, periods, periods_errs = get_periods()

    #dist_pc, dist_pc_errs = get_parallaxes()

    dist_pc, dist_pc_errs, bvcorrection, bverror = read_parallaxes()

    objs, periods, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier = plot_calc_rm(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection, bverror)

    slope_array, offset_array = jackknifing_a(abs_mags, abs_mags_errs, periods)

    jacky = jackgen(slope_array, offset_array, log_periods)

    plot_pl_rm(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, outlier, jacky)

    aliasing_rm(log_periods, abs_mags, abs_mags_errs, slope_array, offset_array)

def pl_gen_odr():
    mags, mags_errs, periods, periods_errs = get_periods()
    #dist_pc, dist_pc_errs = get_parallaxes()
    dist_pc, dist_pc_errs, bvcorrection, bverror = read_parallaxes()
    objs, periods, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, beta, slope_array, offset_array = orthogonal_distance_regression(mags, mags_errs, periods, periods_errs, dist_pc, dist_pc_errs, bvcorrection, bverror)
    plot_pl_odr(objs, abs_mags, abs_mags_errs, log_periods, log_periods_errs, diff_inerrs, beta)
    aliasing_odr(log_periods, abs_mags, log_periods_errs, abs_mags_errs, slope_array, offset_array)

#pl_gen()
#pl_gen_sans_outliers()
#pl_gen_odr()
