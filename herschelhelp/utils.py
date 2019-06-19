# -*- coding: utf-8 -*-

import logging
import os

import numpy as np

import pyvo as vo


    
def clean_table(table):
    """Take a table produced by a VO query and remove all empty columns
    
    Often many columns are empty and make the tables hard to read.
    The function also converts columsn that are objects to strings.
    Object columns prevent writing to fits.
    
    Inputs
    =======
    table
    
    Returns
    =======
    table
    
    """
    table = table.copy()
    for col in table.colnames:
        if np.all(table[col].mask):
            table.remove_column(col)
    for col in table.colnames:
        if table[col].dtype == 'object':
            print(col, table[col].dtype )
        
            table[col] = table[col].astype(str)
    return table