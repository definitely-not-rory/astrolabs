from imports import *
from modes_plot import modes_plot

def overfitted_modes(objs, periods):
    pos=0
    overfitted_modes=[]
    for obj in objs:
        overfitted_modes.append(modes_plot(obj,periods[pos],10,False))
        pos+=1
    plt.figure()
    plt.hist(overfitted_modes)
    plt.xlabel('Fourier coefficient $n$ where $\chi^2<1$ in fitting for $\sum_{n=1}^{10}A_n\sin(n\\frac{\pi}{T}t+\phi_n)$')
    plt.show()
    