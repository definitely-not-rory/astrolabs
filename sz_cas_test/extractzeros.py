import numpy as np
def get_zeros():
    zeros = []
    err_zeros = []

    with open('zeros.txt','r') as readfile:
        for line in readfile:
            line = line.split()
            zeros = np.append(zeros,line[2])
            err_zeros = np.append(err_zeros,line[5])

    return zeros, err_zeros
        
        
