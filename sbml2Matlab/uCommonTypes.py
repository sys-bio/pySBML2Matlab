# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 21:19:31 2020

@author: hsauro
"""

from dataclasses import dataclass

@dataclass
class SpeciesData:
    speciesName : str = ''  # Not yet collected
    speciesId  : str = ''         
    boundaryCondition : bool = True
    initialConcentration : float = 0.0
    isConcentration : bool = True 
    initialAmount : float = 0.0
    isAmount : bool = False
    compartmentId  : str = ''
    compartmentVolume : float = 0.0

class MapElement:
    
    def __init__(self, index, matlabSymbol):
       self.index = index  # Index refers to the matlab indexed array
       self.matlabSymbol = matlabSymbol
    
@dataclass
class NameValue:
	Id    : str = ''
	value : float = 0.0