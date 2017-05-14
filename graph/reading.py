# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 16:17:46 2016

@author: claudio
"""

def read_from_file(name):
    """Function to read from a file called name

    it return a matrix of the various numbers"""
    f= open(name,'r')
    lines=f.readlines()
    lines2=[x.replace("\t"," ") for x in lines]
    lines3=[x.replace("\n","") for x in lines2]
    result=[[float(x) for x in res.split()] for res in lines3]
    return result
    
def obtain_x(name1,column1):
    """Function that exctract a vectors of numbers from a file
    
    it return the x (or y) vector"""
    result1=read_from_file(name1)
    x=[ris1[column1-1] for ris1 in result1]
    return x

    
    