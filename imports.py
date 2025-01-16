import numpy as np
import pandas as pd
import scipy as sp
import astropy as astro
from astropy import time as astrotime
import matplotlib
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
import os
import math
import re
import datetime
import time
import random
import sys
from ast import literal_eval
import warnings
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from scipy.stats import norm 
import statistics 
plt.rcParams.update({'font.size': 14})
warnings.filterwarnings("ignore")

import matplotlib.font_manager as font_manager
plt.rcParams["axes.formatter.use_mathtext"]=True
plt.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=matplotlib.get_data_path() + '/fonts/ttf/cmr10.ttf')
plt.rcParams['font.serif']=cmfont.get_name()
plt.rcParams['mathtext.fontset']='cm'
plt.rcParams['axes.unicode_minus']=False

### Use pre-release version of astroquery

### pip install -U --pre astroquery
