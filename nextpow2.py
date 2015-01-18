# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 16:07:05 2013

@author: KyoungWon
"""

def nextpow2(i):
    n=2
    while n<i:
        n=n*2
    return n
    