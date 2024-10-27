#IMPORTS
import numpy as np
import pandas as pd
import scipy as sp
import astropy as astro
import matplotlib.pyplot as plt
import os
import math

#DATA PROCESSING
def get_data(obj):
    dates =os.listdir(obj+'/') #get list of available nights for given object
    times=[] #storage for data
    mags=[]
    errors=[]
    try: #try/except statement to ignore non-directory files/photometry files
        for date in dates: #iterate through available data
            df = pd.read_csv(obj+'/'+date+'/results.diff', delimiter=' ') #read results file outputted from raw2dif.py
            arrays=df.to_numpy()[:,:3] #remove NaN values (idk why they're there)
            for point in arrays: #iterate through extracted data and store in correct places
                times.append(point[0])
                mags.append(point[1])
                errors.append(point[2])
    except:
        pass #ignore all non-directories
    times-=times[0] #normalise time values to begin at 0
    return times, mags, errors


#INITIAL PLOTS
def plot(objs):
    rows=math.ceil(len(objs)/3) #obtain rows needed to generate 3xhowever many grid
    if rows>1: #Since plots < 3 are stored in a 1D array, need to handle plotting separately
        figure, ax = plt.subplots(rows,3)  #generate subplots
        figure.set_figwidth(15)
        figure.set_figheight(rows*10) #setting figure dimensions

        xcounter=0 #counters for iterations
        ycounter=0

        xaxiscounter=0
        yaxiscounter=0

        for axis in range(rows*3):
            if (xaxiscounter+1)+yaxiscounter*3>len(objs): #check if 'integer' value is within necessary list of objects
                ax[yaxiscounter][xaxiscounter].set_visible(False) #hide axis if plot is unecessary
            if xaxiscounter==2: #resetting x position if at end of the row. increment y position to next row
                xaxiscounter=0
                yaxiscounter+=1
            else:
                xaxiscounter+=1 #if in middle of row, progress to next x position

        for obj in objs:
            times, mags, errors =get_data(obj) #get data from given object
            markers,bars,caps=ax[ycounter][xcounter].errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=2) #generate figure

            [bar.set_alpha(0.5) for bar in bars] #set error bars to translucent
            [cap.set_alpha(0.5) for cap in caps] 

            ax[ycounter][xcounter].set(xlabel='Time (days)',ylabel='Magnitude') #label with axes labels and relevant object name
            ax[ycounter][xcounter].set_title(obj)

            if xcounter==2: #same resetting of x and y positions for row iteration as above
                xcounter=0
                ycounter+=1
            else:
                xcounter+=1
    else:
        figure, ax = plt.subplots(rows,3)  #identical code without y positions for 1 row plots
        figure.set_figwidth(60) 
        
        xcounter=0
        xaxiscounter=0
        for axis in ax:
            if xaxiscounter+1>len(objs):
                ax[xaxiscounter].set_visible(False)
            xaxiscounter+=1

        for obj in objs:
            times, mags, errors =get_data(obj)
            markers,bars,caps=ax[xcounter].errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=3)
            [bar.set_alpha(0.5) for bar in bars]
            [cap.set_alpha(0.5) for cap in caps]
            ax[xcounter].set(xlabel='Time (days)',ylabel='Magnitude')
            #ax[xcounter].set_title(obj)
            xcounter+=1
    
    plt.show()

#CURVE FITTING
def fitting(obj):
    times, mags, errors =get_data(obj) #get data for given object

    initial_values = [max(mags)-(max(mags)+min(mags))/2,0.5,2*np.pi,(max(mags)+min(mags))/2] #setting of trial values for both models
    
    bounds = ([-np.inf,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf])

    def sin_function(t,*params):
        '''
        amplitude = params[0]
        period = params[1]
        phase = params[2]
        displacement = params[3]
        '''
        return params[0]*np.sin(params[1]*t-params[2])+params[3] #defining sine function for fitting
   
    def sawtooth_function(t, *params):
        return params[0]*sp.signal.sawtooth(params[1]*t-params[2])+params[3] #defining sawtooth function for fitting
    

    sawpopt, sawcov = sp.optimize.curve_fit(sawtooth_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6, bounds=bounds) #run fitting for each model
    sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6, bounds=bounds)
    
    smooth_x=np.linspace(times[0], times[-1], 1000) #define x-range for plotting

    plt.errorbar(times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3) #display raw data

    plt.plot(smooth_x,sawtooth_function(smooth_x, *sawpopt),c='r',linestyle='dashed') #plot fitted models
    plt.plot(smooth_x,sin_function(smooth_x, *sinpopt),c='b',linestyle='dashed')

    plt.xlabel('Time (days)') #axes errors
    plt.ylabel('Magnitude')

    sawpopt_errs = np.sqrt(np.diag(sawcov)) #calculate errors
    sinpopt_errs = np.sqrt(np.diag(sincov))

    print(sinpopt[1])
    print(sawpopt[1])

    print('Sawtooth Freq: '+str(sawpopt[1]))
    print('Sin Period: '+str(sinpopt[1])+'\n')

    print('Sawtooth Period: ' +str(2*np.pi/sawpopt[1])+ ' +/- '+str(sawpopt_errs[1]/sawpopt[1]**2)) #print calculated periods with errors
    print('Sin Period: '+str(2*np.pi/sinpopt[1])+' +/- '+str(sinpopt_errs[1]/sinpopt[1]**2))

    plt.show()


objs=['ch_cas','cg_cas','sz_cas','cp_cep','cy_cas']

#plot(objs)

def basic_plot(obj):
    fig, ax=plt.subplots()
    times, mags, errors =get_data(obj)
    markers,bars,caps=ax.errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=3)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    ax.set(xlabel='Time (days)',ylabel='Magnitude')
    fig.set_figheight(5)
    fig.set_figwidth(7.5)
    plt.title(obj)
    plt.show()



#basic_plot('cg_cas')

def fittingtesting():
    times, mags, errors = get_data('cg_cas')
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
    def chi_squared(model_params, model, x_data, y_data, y_err):
        chi_sq = np.sum(((y_data - model(x_data, model_params))/y_err)**2)
        return chi_sq
    

    print(chi_squared(fitted_params,linear,times,mags,errors))
    

def aphot_cg_cas_testing():
    times, mags, errors = get_data('cg_cas')
    print(str(mags[0])+' +/- '+str(errors[0]))
    print(str(mags[1])+' +/- '+str(errors[1]))


fittingtesting()

'''

chi squared values

errors from chi+1

check data errors

test at single points rather than every exposure

average and sd of errors

include errors in ZP

do github

CG CAS CHECK:

ad101.fits = 11.758 +/- 0.019
ad102.fits = 11.725 +/- 0.018


'''