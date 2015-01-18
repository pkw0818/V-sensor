# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 12:52:26 2013

@author: KyoungWon
"""

import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt, numpy as np, numpy.random, scipy

data = genfromtxt('C:\Users\KyoungWon\Google Drive\V sensor\dIoverdE.csv', delimiter = ',') 
#dE=hist(data[:,0])
#dI=hist(data[:,1])
#dIvalue=hist(data[:,0])
dE=dE[1]
dI=dI[1]

scatter(dE,dI,c=data[:,0])