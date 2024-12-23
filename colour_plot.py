from imports import *
from fourier import fourier_fitting

def colour_plot(obj,period,n):
    n1s=range(1,n+1)
    n2s=range(1,n+1)
    chis=np.zeros((n,n))

    for n1 in n1s:
        for n2 in n2s:
            chis[n1-1][n2-1]=fourier_fitting(obj,period,n1,n2,False,False,20)[2]

    plt.figure()
    im=plt.imshow(chis,cmap='plasma',extent=(1,n,1,n),origin='lower',norm='log')
    plt.xlabel('$n_1$ coefficient in 2 mode Fourier function')
    plt.ylabel('$n_2$ coefficient in 2 mode Fourier function')
    plt.colorbar(im,label='Reduced $\chi^2$ value of $n_1$, $n_2$ 2 mode Fourier function')
    plt.show()
