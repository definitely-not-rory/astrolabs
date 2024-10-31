from imports import *
from astrolabs import *

def plot_grid(objs):
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
            initial_values = [max(mags)-(max(mags)+min(mags))/2,0.5,2*np.pi,(max(mags)+min(mags))/2] #setting of trial values for both models
            def sin_function(t,*params):
                '''
                amplitude = params[0]
                frequency = params[1]
                phase = params[2]
                displacement = params[3]
                '''
                return params[0]*np.sin(params[1]*t-params[2])+params[3] #defining sine function for fitting
            
            sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6)
            smooth_x=np.linspace(times[0], times[-1], 1000) #define x-range for plotting
            
            markers,bars,caps=ax[ycounter][xcounter].errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=2) #generate figure
            
            ax[ycounter][xcounter].plot(smooth_x,sin_function(smooth_x, *sinpopt),c='b',linestyle='dashed')
            
            #[bar.set_alpha(0.5) for bar in bars] #set error bars to translucent
            #[cap.set_alpha(0.5) for cap in caps] 

            ax[ycounter][xcounter].set(xlabel='Time (days)',ylabel='Magnitude') #label with axes labels and relevant object name
            ax[ycounter][xcounter].set_title(obj)

            sinpopt_errs = np.sqrt(np.diag(sincov))

            def chi_squared(model_params, model, x_data, y_data, y_err):
                return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
            sin_chi_val=chi_squared(sinpopt, sin_function, times, mags, errors)
            reduced_sin_chi=sin_chi_val/len(times)

            print('\n~~~ '+obj+' Sinusoidal Model ~~~\nSin Frequency: '+str(sinpopt[1]))
            print('Sin Period: '+str(2*np.pi/sinpopt[1])+' +/- '+str(sinpopt_errs[1]/sinpopt[1]**2))
            print('Sinusoidal Reduced Chi Squared: '+str(reduced_sin_chi)+'\n')

            if xcounter==2: #same resetting of x and y positions for row iteration as above
                xcounter=0
                ycounter+=1
            else:
                xcounter+=1
        figure.tight_layout()
    else:
        figure, ax = plt.subplots(rows,3)  #identical code without y positions for 1 row plots
        figure.set_figwidth(60) 
        figure.tight_layout()

        xcounter=0
        xaxiscounter=0
        for axis in ax:
            if xaxiscounter+1>len(objs):
                ax[xaxiscounter].set_visible(False)
            xaxiscounter+=1

        for obj in objs:
            times, mags, errors =get_data(obj)
            markers,bars,caps=ax[xcounter].errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=3)
            #[bar.set_alpha(0.5) for bar in bars]
            #[cap.set_alpha(0.5) for cap in caps]

            initial_values = [max(mags)-(max(mags)+min(mags))/2,0.5,2*np.pi,(max(mags)+min(mags))/2] #setting of trial values for both models
            def sin_function(t,*params):
                '''
                amplitude = params[0]
                frequency = params[1]
                phase = params[2]
                displacement = params[3]
                '''
                return params[0]*np.sin(params[1]*t-params[2])+params[3] #defining sine function for fitting
            
            sinpopt, sincov = sp.optimize.curve_fit(sin_function,times,mags,sigma=errors,absolute_sigma=True,p0=initial_values,check_finite=True, maxfev=10**6)
            smooth_x=np.linspace(times[0], times[-1], 1000) #define x-range for plotting

            ax[xcounter].set(xlabel='Time (days)',ylabel='Magnitude')
            ax[xcounter].plot(smooth_x,sin_function(smooth_x, *sinpopt),c='b',linestyle='dashed')

            sinpopt_errs = np.sqrt(np.diag(sincov))

            def chi_squared(model_params, model, x_data, y_data, y_err):
                return np.sum(((y_data - model(x_data, *model_params))/y_err)**2)
    
            sin_chi_val=chi_squared(sinpopt, sin_function, times, mags, errors)
            reduced_sin_chi=sin_chi_val/len(times)

            print('\n~~~ '+obj+' Sinusoidal Model ~~~\nSin Frequency: '+str(sinpopt[1]))
            print('Sin Period: '+str(2*np.pi/sinpopt[1])+' +/- '+str(sinpopt_errs[1]/sinpopt[1]**2))
            print('Sinusoidal Reduced Chi Squared: '+str(reduced_sin_chi)+'\n')

            ax[xcounter].set_title(obj)
            xcounter+=1
    
    plt.show()

objs=['cg_cas','ch_cas','cp_cep','cy_cas','sz_cas','kx_cyg','v396_cyg','v532_cyg','v538_cyg','v609_cyg','vx_cyg']

plot_grid(objs)

