# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 13:00:30 2013

@author: KyoungWon
"""
#
def clearall():
    """clear all globals"""
    for uniquevar in [var for var in globals().copy() if var[0] != "_" and var != 'clearall']:
        del globals()[uniquevar]