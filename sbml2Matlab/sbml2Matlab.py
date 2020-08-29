# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 17:43:35 2020

@author: hsauro
"""

import simplesbml
import warnings

# try to import tesbml or libsbml
# if both of these fail, libsbml cannot be imported - cannot continue
try:
    import tesbml as libsbml   
except ImportError:
    import libsbml

from math import isnan
from re import sub
import os

# Version info is in __init__.py


class sbml2Matlab(object):
      pass
    