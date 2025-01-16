from imports import *
from new_fourier import new_fourier

def gen_pmf(objs,periods):
    f=open('periodmagdata.txt','w')
    f.write('')
    f.close()

    f=open('fitting_params.txt','w')
    f.write('obj lit_period max_modes max_period max_chi n1 n2 bimodal_period bimodal_errs bimodal_chi\n')
    f.close()

    f=open('plot_data.txt','w')
    f.write('')
    f.close()
    pos=0
    for obj in objs:
        correct=False
        while correct==False:
            m_0, m_0_err,p,p_err,fitting_params=new_fourier(obj,periods[pos])
            correct=bool(input('Good Fit?: '))
        f = open("periodmagdata.txt", "a")
        f.write(str(m_0)+'x'+str(m_0_err)+'x'+str(p)+'x'+str(p_err)+'\n')
        f.close()

        f=open('fitting_params.txt','a')
        for i in fitting_params:
            f.write(str(i)+' ')
        f.write('\n')
        f.close()
        pos+=1
