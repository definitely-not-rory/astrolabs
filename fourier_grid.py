from imports import *
from data_handling import get_data
from fourier import fourier_fitting

def plot_fourier_grid(objs, periods,folded):
    f=open('periodmagdata.txt','w')
    f.write('')
    f.close()
    
    df = pd.read_csv("mcmaster.txt",delimiter=",")
    averagemags = df["mean_mag"].to_numpy()

    rows=math.ceil(len(objs)/3) #obtain rows needed to generate 3xhowever many grid
    if rows>1: #Since plots < 3 are stored in a 1D array, need to handle plotting separately
        figure, ax = plt.subplots(rows,3)  #generate subplots
        figure.set_figwidth(15)
        figure.set_figheight(rows*10) #setting figure dimensions
        figure.tight_layout()

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

        obj_counter=0
        for obj in objs:

            times, mags, errors,days =get_data(obj) #get data from given object

            period=periods[obj_counter]

            #mean_mag=np.mean(mags)
            mean_mag = averagemags[obj_counter]
            #mean_mag_error=np.std(mags)/(np.sqrt(len(mags)))
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
            p_lo=period/1.2
            p_hi=period*1.2
            phi1_lo=0
            phi1_hi=2*np.pi
            phi2_lo=0
            phi2_hi=2*np.pi
            disp_lo=mean_mag-1
            disp_hi=mean_mag+1

            fourier_bounds=([disp_lo,amp1_lo,p_lo,phi1_lo,amp2_lo,phi2_lo],[disp_hi,amp1_hi,p_hi,phi1_hi,amp2_hi,phi2_hi])
            
            jackknifed_mean_mags = []
            for i in range(len(times)):
                jackknifed_times = np.delete(times,i)
                jackknifed_mags = np.delete(mags,i)
                jackknifed_errors = np.delete(errors,i)
                popt,cov=sp.optimize.curve_fit(fourier_function,jackknifed_times,jackknifed_mags,sigma=jackknifed_errors,p0=fourier_values,bounds=fourier_bounds,check_finite=True,maxfev=10**6)
                jackknifed_mean_mags = np.append(jackknifed_mean_mags,popt[0])
                
            popt,cov=sp.optimize.curve_fit(fourier_function,times,mags,sigma=errors,p0=fourier_values,bounds=fourier_bounds,check_finite=True,maxfev=10**6)
            
            mean_mag = popt[0]
            mean_mag_error = np.std(jackknifed_mean_mags)

            smooth_x=np.linspace(times[0], times[-1], 1000)

            fitted_period,fitted_error,reduced_chi,mean_mag,mean_mag_error=fourier_fitting(obj,period,2,5,False,False,20)

            f = open("periodmagdata.txt", "a")
            f.write(str(mean_mag)+' '+str(mean_mag_error)+' '+str(fitted_period)+' '+str(fitted_error)+'\n')
            f.close()

            details=obj+f'\nFitted Period: {fitted_period:.2f}+/-{fitted_error:.1g}\nReduced Chi Squared: {reduced_chi:.2f}'

            if folded==False:
                markers,bars,caps=ax[ycounter][xcounter].errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=2)

                ax[ycounter][xcounter].plot(smooth_x,fourier_function(smooth_x, *popt),c='b',linestyle='dashed')

                ax[ycounter][xcounter].set(xlabel='Time (days)',ylabel='Magnitude') #label with axes labels and relevant object name
                ax[ycounter][xcounter].invert_yaxis()
                ax[ycounter][xcounter].text(min(ax[ycounter][xcounter].get_xlim())+0.25,max(ax[ycounter][xcounter].get_ylim())+0.05,details)
            else:
                folded_times=times
                for i in range(len(folded_times)):
                    folded_times[i]=folded_times[i]%fitted_period

                folded_fit_times=np.linspace(0,fitted_period,1000)
                ax[ycounter][xcounter].errorbar(folded_times,mags,yerr=errors,marker='x',linestyle='None',c='k',capsize=3)
                ax[ycounter][xcounter].plot(folded_fit_times,fourier_function(folded_fit_times,*popt),c='r')
                ax[ycounter][xcounter].invert_yaxis()
                ax[ycounter][xcounter].text(min(ax[ycounter][xcounter].get_xlim()),max(ax[ycounter][xcounter].get_ylim()),details)

            if xcounter==2: #same resetting of x and y positions for row iteration as above
                xcounter=0
                ycounter+=1
            else:
                xcounter+=1
            
            obj_counter+=1
    else:
        figure, ax = plt.subplots(rows,3)  #identical code without y positions for 1 row plots
        figure.set_figwidth(60) 
        figure.tight_layout()

        xcounter=0
        xaxiscounter=0
        obj_counter=0

        for axis in ax:
            if xaxiscounter+1>len(objs):
                ax[xaxiscounter].set_visible(False)
            xaxiscounter+=1

        for obj in objs:
            times, mags, errors,days =get_data(obj) #get data from given object

            period=periods[obj_counter]

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
            p_lo=period/1.2
            p_hi=period*1.2
            phi1_lo=0
            phi1_hi=2*np.pi
            phi2_lo=0
            phi2_hi=2*np.pi
            disp_lo=np.mean(mags)-1
            disp_hi=np.mean(mags)+1

            fourier_bounds=([disp_lo,amp1_lo,p_lo,phi1_lo,amp2_lo,phi2_lo],[disp_hi,amp1_hi,p_hi,phi1_hi,amp2_hi,phi2_hi])

            popt,cov=sp.optimize.curve_fit(fourier_function,times,mags,sigma=errors,p0=fourier_values,bounds=fourier_bounds,check_finite=True,maxfev=10**6)

            smooth_x=np.linspace(times[0], times[-1], 1000)

            markers,bars,caps=ax[xcounter].errorbar(times,mags,errors,fmt='o',c='r', marker='x',ecolor='k',capsize=2)

            ax[xcounter].plot(smooth_x,fourier_function(smooth_x, *popt),c='b',linestyle='dashed')

            ax[xcounter].set(xlabel='Time (days)',ylabel='Magnitude') #label with axes labels and relevant object name

            xcounter+=1
            obj_counter+=1
    plt.show()

    
