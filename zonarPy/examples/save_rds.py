# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 09:45:47 2020

@author: sveng
"""
from rpy2 import robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()

def save_rds(df, filename):
    r_data = pandas2ri.py2ri(df)
    robjects.r.assign("my_df", r_data)
    robjects.r("saveRDS(my_df, file='{}', comp='gzip')".format(filename))