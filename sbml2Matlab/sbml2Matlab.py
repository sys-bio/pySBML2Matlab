# -*- coding: utf-8 -*-
"""

@author: hsauro
"""

# sbml2matalb convertter
#  Does not support events, local parameters or algebraic rule

# Questions: 
   # Are initialization expression initialized?
   
import simplesbml
import warnings
from   uSubsitutionCode import *
from   uCommonTypes import *

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


def _isSBMLModel(obj):
  """
  Tests if object is a libsbml model
  """
  cls_stg = str(type(obj))
  if ('Model' in cls_stg) and ('lib' in cls_stg):
    return True
  else:
    return False

def _checkSBMLDocument(document): 
  if (document.getNumErrors() > 0):
    raise ValueError("Errors in SBML document")

    
class ReactionInfo:
    
       def __init__(self, sbmlModel : simplesbml.SbmlModel, reactionIndex : int):

        self.Id : str
        self.name : str  # Not yet used
        self.rateLaw : str
    
        self.reactants  = []  # NameValue
        self.products   = []  # NameValue
        self.parameters = []  # parameters
      
        self.Id = sbmlModel.getNthReactionId (reactionIndex)    
        self.rateLaw = sbmlModel.getRateLaw (reactionIndex)
    
        numOfReactants = sbmlModel.getNumReactants (reactionIndex)
        for i in range (numOfReactants):
            reactant = NameValue()
            reactant.Id = sbmlModel.getReactant (reactionIndex, i)
            reactant.value = sbmlModel.getReactantStoichiometry (reactionIndex, i)
            self.reactants.append (reactant)
    
        numOfProducts = sbmlModel.getNumProducts (reactionIndex);
        for i  in range (numOfProducts):
            reactant = NameValue()
            reactant.Id = sbmlModel.getProduct (reactionIndex, i)
            reactant.value = sbmlModel.getProductStoichiometry (reactionIndex, i)
            self.products.append (reactant)
    
        # Not yet implemented
        # numParameters = getNumLocalParameters (reactionIndex)
        # #     for (int i=0; i<numParameters; i++) {
        # #         TNameValue parameter;
        # #         getNthLocalParametesperId (reactionIndex, i, &cId);
        # #         parameter.name = cId;
        # #         getNthLocalParameterValue (reactionIndex, i, &value);
        # #         parameter.value = value;
        # #         parameters.push_back(parameter);
        # #     }
        # # }
        
        
class Sbml2Matlab(object):
    
    def __init__(self, sbmlStr=None, sbmlFile=None):
        
        self.model = None
        if sbmlStr != None:
           self.model = simplesbml.SbmlModel(sbmlStr=sbmlStr)
        if sbmlFile != None:
           self.model = simplesbml.SbmlModel (sbmlFile=sbmlFile)   
    
        self.speciesData : SpeciesData
        self.reactions : ReactionInfo
        self.parameterMap = {}
        self.speciesMap = {}
        
    def convert (self):
        self._readSpecies()
        self._readParameters()
        self._readReactions()
        
        # Avoid underscores
        if self.model.getModelId() == '__main':
            self.modelId = 'main'
        else:
            self.modelId = self.model.getModelId()
            
        result  = self.printHeader()
        result += self.printSpeciesOverview()
        result += self.printOutCompartments()
        result += self.printOutGlobalParameters() 
        result += self.printOutEventFlags()
        result += self.printOutBoundarySpecies()
        result += self.printInitialConditions()
        result += "\nelse\n" # End of if (nargin == 0) block
        result += self.printOutRules()
        result += self.printOutEvents()          
        result += self.printRatesOfChange()           
        result += self.printOutReactionScheme()       
        result += self.printSupportedFunctions()
        
        return result
       
    # Method that creates one Species object
    # for each species encapsulating all available information. This includes
    # the following
    # - species name          string
    # - species id            string
    # - boundaryCondition     bool
    # - initialConcentration  double
    # - initialAmount         double
    # - compartment           string
    # - compartmentVolume     double
    def _readSpecies(self): 

        numTotalSpecies = self.model.getNumFloatingSpecies() + self.model.getNumBoundarySpecies()

        for i in range (self.model.getNumFloatingSpecies()): 
            matlabSymbol = 'x(' + str (i+1) + ')'
            self.speciesMap[self.model.getNthFloatingSpeciesId (i)] = MapElement (i, matlabSymbol)

        self.speciesData = []

        for i in range (self.model.getNumFloatingSpecies()):
            spElement = SpeciesData()
            self.speciesData.append (spElement)
            
            spElement.speciesId = self.model.getNthFloatingSpeciesId(i)

            spElement.compartmentId = self.model.getCompartmentIdSpeciesIsIn (spElement.speciesId)  
            spElement.compartmentVolume =  self.model.getCompartmentVolume (spElement.compartmentId)
            
            if self.model.isConcentration (spElement.speciesId):
               spElement.initialConcentration = self.model.getSpeciesInitialConcentration(spElement.speciesId)
               spElement.isConcentration = True
               spElement.initialAmount = spElement.initialConcentration*spElement.compartmentVolume
               
            if self.model.isAmount (spElement.speciesId):
               spElement.initialAmount = self.model.getSpeciesInitialAmount(spElement.speciesId)
               spElement.initialConcentration = spElement.initialAmount/spElement.compartmentVolume
               spElemet.isAmount = True

            # To be done 
            #getNthFloatingSpeciesName(i, &cstr);
            #sp_list[i].name = cstr;
            
        for i in range (self.model.getNumBoundarySpecies()):

            spElement = SpeciesData()
            self.speciesData.append (spElement)             
            
            spElement.speciesId = self.model.getNthBoundarySpeciesId (i)

            spElement.compartmentId = self.model.getCompartmentIdSpeciesIsIn (spElement.speciesId)  
            spElement.compartmentVolume =  self.model.getCompartmentVolume (spElement.compartmentId)

            if self.model.isConcentration (spElement.speciesId):
               spElement.initialConcentration = self.model.getSpeciesInitialConcentration(spElement.speciesId)
               spElement.isConcentration = True
               spElement.initialAmount = spElement.initialConcentration*spElement.compartmentVolume
               
            if self.model.isAmount (spElement.speciesId):
               spElement.initialAmount = self.model.getSpeciesInitialAmount()
               spElement.initialConcentration = spElement.initialAmount/spElement.compartmentVolume
               spElemet.isAmount = True

    def _readParameters (self):

        # Read the parameters (including boundary species) and map them to indices
        for i in range (self.model.getNumParameters()): 
            matlabSymbol = 'gp_' + str (i+1)
            self.parameterMap[self.model.getParameterId (i)] = MapElement (i, matlabSymbol)
        for i in range (self.model.getNumBoundarySpecies()):
            matlabSymbol = 'gp_' + str (i+self.model.getNumParameters()+1)
            self.parameterMap[self.model.getNthBoundarySpeciesId (i)] = MapElement (i+self.model.getNumParameters(), matlabSymbol)

    # Collection information on each reaction
    def _readReactions(self):
        self.reactions = []
        numReactions = self.model.getNumReactions()
        for i in range (numReactions):
            self.reactions.append(ReactionInfo(self.model, i))
    
    # prints the header information on how to use the matlab file
    def printHeader(self):
        result = '' 
        result += "function xdot = " + self.modelId + "(time,x)\n"
        result +=  "%  synopsis:\n"
        result +=  "%     xdot = " + self.modelId + " (time, x)\n"
        result +=  "%     x0 = " + self.modelId +"\n"
        result +=  "%\n"
        result +=  "%  " + self.modelId  + " can be used with the Matlab ode functions as follows:\n%\n"
        result +=  "%  x0 = " + self.modelId + ";\n"
        result +=  "%  [t,x] = ode23s(@" + self.modelId + ", [0 100], " + self.modelId + ");\n"
        result +=  "%  plot (t,x);\n"
        result +=  "%\n"
        result +=  "%  where 100 is the end time\n" + "%\n"
        result +=  "%  When " + self.modelId + " is used without any arguments it returns a vector of\n"
        result +=  "%  the initial concentrations of the the floating species.\n"
        result +=  "%  Otherwise " + self.modelId + " should be called with two arguments:\n"
        result +=  "%  time and x.\n"
        result +=  "%     time is the current time.\n"
        result +=  "%     x is the vector of the concentrations of the " + str (self.model.getNumFloatingSpecies ()) + " floating species.\n"
        result +=  "%  When these parameters are supplied, " + self.modelId + " returns a vector of\n"
        result +=  "%  the rates of change of the concentrations of the " + str (self.model.getNumFloatingSpecies ()) + " floating species.\n" + "%\n"

        result +=  "%  NOTE for compartmental models the matlab translator\n"
        result +=  "%  generates code that when simulated in matlab, produces\n"
        result +=  "%  results which have the units of species amounts. Users\n"
        result +=  "%  should divide the results for each species with the volume of the\n"
        result +=  "%  compartment it resides in, in order to obtain concentrations.\n"
        result +=  "%  \n"
        
        result +=  "%  The following table shows the mapping between the vector\n"
        result +=  "%  index of a floating species and the species name.\n\n" 

        result +=  "%  Indx      Name\n"
        return result
    
    def printSpeciesOverview (self):
        result = ''
        for i in range (self.model.getNumFloatingSpecies()):      
            result +=  "%  x(" + str(i+1) +  ")        " + self.speciesData[i].speciesId + '\n'
        return result
     
    # prints out the compartment information
    def printOutCompartments(self):

        result = "\n% List of Compartments\n"       

        compartmentList = self.model.getListOfCompartmentIds ()
        for i in range (self.model.getNumCompartments ()):
            result += "vol__" + compartmentList[i] + " = " + str (self.model.getCompartmentVolume(i))  + ';\n'
        #   + ";\t\t%"  + compartmentList[i].name + "\n"
     
        return result

    # prints out the list of global parameters
    def printOutGlobalParameters(self):
    
       result = ''      

       if self.model.getNumParameters() > 0:
            result += "\n% Global Parameters\n"
       for i in range (self.model.getNumParameters()):
           result +=  "gp_" + str (i+1) + " = " + str (self.model.getParameterValue (i)) + ";\t\t% " \
              + self.model.getParameterId (i) + '\n'

       return result


    # prints out the boundary species
    def printOutBoundarySpecies(self):
        
        result = ''
        if self.model.getNumBoundarySpecies () > 0:
            result += "\n% Boundary Conditions\n"

        for i in range (self.model.getNumBoundarySpecies ()):
            index = i+self.model.getNumFloatingSpecies()
            speciesId = self.speciesData[index].speciesId
            isAmount = self.model.isAmount (speciesId)

            result +=  "rInfo.g_p" + str (self.model.getNumParameters() + i+1) + " = ";             

            if (isAmount == True):
               value =  self.model.getSpeciesInitialAmount (speciesId)                
            else:
               value = self.model.getSpeciesInitialConcentration(speciesId)                

            result += str (value) + ";\t\t% " + speciesId + '\n' 
            ##<< (isAmount ? " [Amount]" : "[Concentration]")  << endl;
 
        return result

    # prints out the initial conditions and reaction info
    def printInitialConditions(self):

        result = "\nif (nargin == 0)" + "\n\n"
        result += "   % set initial conditions\n"

        initCondIndex = 1; # index of current initial condition
        for i in range (self.model.getNumFloatingSpecies()):

            floatingSpeciesName = self.speciesData[i].speciesId
            bnd_data = ''
            if self.speciesData[i].isAmount == False:
               value = self.speciesData[i].initialAmount
               bnd_data = "Concentration"
               strValue = str (value)
            else:
              value = self.speciesData[i].initialConcentration
              bnd_data = "Amount converted to concentration"
              strValue = str (value)
              strValue = strValue + "/vol__" + self.speciesData[i].compartmentId

            result +=  "   xdot(" + str (i+1) + ") = " + strValue  + ";\t% " + floatingSpeciesName \
                  + " = " + self.speciesData[i].speciesName + bnd_data + '\n'
            initCondIndex += 1

        return result

    # ---------------------------------------------------------
    # Print out the event flag 
    def printOutEventFlags (self):
        result = ''
        if self.model.getNumEvents():
           result += '\n% Event flags\n' 
           for i in range (self.model.getNumEvents()):
               result += 'T_' + self.model.getEventId (i) + ' = true;\n'
                  
        return result           
            
        # ---------------------------------------------------------
        # printing out boundary species
        if self.model.getNumBoundarySpecies() > 0:
            result += "\n   rInfo.boundarySpecies = {" + " % Each row: [Species Id, Initial Value, isAmount (1 for amount, 0 for concentration)]\n"
            for i in range (self.model.getNumBoundarySpecies()):
                index = i + self.model.getNumFloatingSpecies()
                isAmount = self.speciesData[index].isAmount
                speciesId = self.speciesData[index].speciesId
    
                result +=  "      '" + speciesId + "', ";
    
                if (isAmount == True):
                   value =  self.speciesData[index].initialAmount
                   valAmount = 1;
                else:
                   value = self.speciesData[index].initialConcentration
                   valAmount = 0;
    
                result += str (value) + ", " + str (valAmount) + '\n'
            result += "   };\n"

        # ---------------------------------------------------------
        # Print out assignment rule information
        if self.model.getNumRules() > 0:
            result += "\n   rInfo.assignmentRules = { \t % List of variables involved in a assignment rule\n" 
            for i in range (self.model.getNumRules()): #adding initial condition for rate rule on a non-species
                if self.model.isRuleType_Assignment(i):
                   variable = self.model.getRuleId (i)
                   equation = self.model.getRuleRightSide (i) 
                   variable = subsituteConstants(self, variable, '', False)
                   equation = subsituteConstants(self, equation, '', False)
                   result += "      " + variable + ' = ' + equation + ";\n"
               
        # ---------------------------------------------------------
        # Print out rate rule information
        if self.model.getNumRules():
            result += "\n   rInfo.rateRules = { \t % List of variables involved in a rate rule\n" 
            for i in range (self.model.getNumRules()): #adding initial condition for rate rule on a non-species
                variable = self.model.getRuleId (i)
                equation = self.model.getRuleRightSide (i) 
                if self.model.isRuleType_Rate(i): # if rate rule, then promote parameter into ode (X)
                   variable = self.model.getRuleId (i)
                   equation = self.model.getRuleRightSide (i) 
    
                   variable = subsituteConstants(self, variable, '', False)
                   equation = subsituteConstants(self, equation, '', False)
                   result += "      " + variable
                   result += " = " + equation + ";\n"
            result += "   };\n" 

        return result
        
    def printOutRules(self):
        result = ''
        if self.model.getNumRules() > 0:
            result += '   % Calculate any assignment rules\n'
            for i in range (self.model.getNumRules()):
                if self.model.isRuleType_Assignment (i):
                   variable = self.model.getRuleId (i)
                   equation = self.model.getRuleRightSide (i) 
    
                   variable = subsituteConstants(self, variable, '', False)
                   equation = subsituteConstants(self, equation, '', False)
                   result += "   " + variable + ' = ' + equation + ";\n"
        return result
        
    def printOutEvents(self):
        result = ''
        if self.model.getNumEvents() > 0:
           for i in range (self.model.getNumEvents ()):
               rightside = subsituteConstants(self, self.model.getEventTrigger (i), '', False) 
               triggerId = 'trigger' + self.model.getEventId(i)
               result += '   ' + triggerId + ' = ' + rightside + ';\n'
               result += '   if T_' + self.model.getEventId(i) + ' && ' + triggerId + '\n'
               for j in range (self.model.getNumEventAssignments (i)):
                   result += '      '
                   result += subsituteConstants(self, self.model.getEventVariable (i, j), '', False) + ' = '
                   result += subsituteConstants(self, self.model.getEventAssignment (i, j), '', False) + ';\n'
                   result += '      T_' + self.model.getEventId(i) + '= false;\n'
               result += '   end;\n'
        return result
            
    def printRatesOfChange(self):
        
        result = "\n   % Calculate rates of change\n"

        for i in range (self.model.getNumReactions()):
            kineticLaw = self.model.getRateLaw (i)
            reactionId = 'X'#_currentModel->reactions[i].id;

            result += "   R" + str (i) +" = " + (subsituteConstants (self, kineticLaw, reactionId)) + ";\n"     
        return result
            
    ###
    def printOutReactionScheme (self):
        result = "\n   xdot = [\n" 
        for i in range (self.model.getNumRules()):
           if self.model.isRuleType_Rate (i):
              variable = self.model.getRuleId (i)
              equation = self.model.getRuleRightSide (i) 

              variable = subsituteConstants(self, variable, '', False)
              equation = subsituteConstants(self, equation, '', False)
              result += "      " + equation + "\t% From rate rule\n"        
            
        xdotIndex = 1
        for i in range (self.model.getNumFloatingSpecies()):
            eqn = "     ";

            floatingSpeciesName = self.speciesData[i].speciesId

            for j in range (self.model.getNumReactions()):
                numProducts = len (self.reactions[j].products)

                for k1 in range (numProducts):
                    productName = self.reactions[j].products[k1].Id

                    if (floatingSpeciesName == productName):
                       productStoichiometry = self.reactions[j].products[k1].value

                       if (productStoichiometry != 1):
                          stoich = str (productStoichiometry) + '*'
                       else:
                          stoich = ""
                          eqn = eqn + " + " + stoich + "R" + str (j)

                numReactants = len (self.reactions[j].reactants)
                for k1 in range (numReactants):
                    reactantName = self.reactions[j].reactants[k1].Id
                    if (floatingSpeciesName == reactantName):
                       reactantStoichiometry = self.reactions[j].reactants[k1].value
                       if (reactantStoichiometry != 1):
                          stoich = str (reactantStoichiometry) + "*"
                       else:
                          stoich = ""
                          eqn = eqn + " - " + stoich + "R" + str (j)

            xdotIndex += 1
            result += eqn + '\n'
            
        result += "   ];\n" 
        result += "end;\n\n"
        return result
        
    def printSupportedFunctions (self):
        result = ''
        return result
        
import tellurium as te

# r = te.loada("""
     
# x1' = k1*x1

# k10 := 45+k11
# k60 := k10*2

# J1: $Xo -> S1; k10*Xo - k11*S1; 
# J2: S1 -> S2; k20*S1 - k21*S2; 
# J3: S2 -> S3; k30*S2 - k31*S3; 
# J4: S3 -> S4; k40*S3 - k41*S4; 
# J5: S4 -> S5; k50*S4 - k51*S5; 
# J6: S5 -> $X1; k60*S5 - k61*X1; 

# k1 = 0.1; x1 = 10
# k11 = 0.69
# k20 = 1.03;  k21 = 0.13
# k30 = 1.89;  k31 = 0.10
# k40 = 4.96;  k41 = 0.61
# k50 = 2.88;  k51 = 0.48
# k61 = 0.83
# Xo = 6.00;   X1 = 5.00
# S1 = 1.1; S2 = 2.2; S3 = 3.3; S4 = 4.4; 
# S5 = 5.5; 
# """)

# r = te.loada("""
# S0 + S3 -> S2; k0*S0*S3;
# S3 + S2 -> S0; k1*S3*S2;
# S5 -> S2 + S4; k2*S5;
# S0 + S1 -> S3; k3*S0*S1;
# S5 -> S0 + S4; k4*S5;
# S0 -> S5; k5*S0;
# S1 + S1 -> S5; k6*S1*S1;
# S3 + S5 -> S1; k7*S3*S5;
# S1 -> $S4 + S4; k8*S1;

# S0 = 1; S1 = 1; S2 = 1; S3 = 1; S4 = 2; S5 = 1;
# k0 = 0.30242636485498864
# k1 = 0.4546338198517117
# k2 = 0.14781881391817286
# k3 = 0.8571005959150634
# k4 = 0.45921360072046713
# k5 = 0.7472466467466189
# k6 = 0.43113361966825015
# k7 = 0.38651480156573215
# k8 = 0.5419808745701772

# S4 = 5
# """)

# m = Sbml2Matlab(r.getSBML())
# print (m.convert())


# r = te.loada("""  
#    k1 := 2*sin (time*5) + 2
#    S1 -> S2; k1*S1;
   
#    S1 = 5; S2 = 0; 
# """)

# m = r.simulate (0, 10, 100, ['time', 'S1', 'k1'])
# r.plot()

# m = Sbml2Matlab(r.getSBML())
# print (m.convert())


r = te.loada("""  
   S1 -> S2; k1*S1;
   
   at time > 20: k1 = 0.8
   at time > 10: k1 = 0.5

   S1 = 5; S2 = 0; k1 = 0.1 
""")

m = Sbml2Matlab(r.getSBML())
print (m.convert())

#m = r.simulate (0, 30, 100, ['time', 'S1', 'k1'])
#r.plot()

