from imports import *
from astrolabs import get_data

def fourier_fitting(obj):
    mags,times,errors=get_data(obj)
    