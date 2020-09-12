# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 17:21:02 2020

@author: hsauro
"""

import lexer

# if expression is parenthesized, drop parentheses
def stringInside(astr : str):
    if (astr == "" or len (astr) < 2):
       return astr
    return astr[1:len(astr)-1]

    
# deal with all strings, which could be: 
# - global parameter (under which we also list boundary species)
# - floating species
# - compartment volumes
# - local parameters
# - function names
# TODO: add flux names!!!
def ReplaceStringToken(currentToken : str, reactionId : str, divideVolumes : bool = True):

    replaceString = ''
    innerString = stringInside(currentToken)        
    localParameterId = reactionId + "_" + innerString;

    if ( _currentModel->globalParamIndexList.find ( innerString ) != _currentModel->globalParamIndexList.end() )
        {
            bool isBoundarySpecies = False

            for ib in range (sbmlModel.getNumBoundarySpecies()):

                if (innerString == _currentModel->sp_list[ib + _currentModel->numFloatingSpecies].id) 
                {
                    if (divideVolumes):
                        replaceStream += "("

                    if (_bInlineMode)
                       replaceStream += _currentModel->globalParametersList[innerString]
                    else:
                        replaceStream += "rInfo.g_p" << _currentModel->globalParamIndexList[innerString]

                    if (divideVolumes):
                    {
                        if (_currentModel->compartmentsList[_currentModel->sp_list[ib + _currentModel->numFloatingSpecies].compartment] != 1.0)
                            replaceStream =+ "/vol__" << _currentModel->sp_list[ib + _currentModel->numFloatingSpecies].compartment;

                        replaceStream += ")"
                    }

                    isBoundarySpecies = True
                    break;
                }
    

            if isBoundarySpecies == False):

                if (_bInlineMode):
                    replaceStream << _currentModel->globalParametersList[innerString]
                else:
                    replaceStream << "rInfo.g_p" << _currentModel->globalParamIndexList[innerString]

            
        }
        else if ( _currentModel->parameterMapList.find ( localParameterId ) != _currentModel->parameterMapList.end() )
        {
            replaceStream << _currentModel->parameterMapList[localParameterId];
        }
        else if ( _currentModel->compartmentsList.find ( innerString ) != _currentModel->compartmentsList.end() )
        {
            replaceStream << "vol__" << innerString;
        }
        else
        {            
            match = False
            for (int isp=0; isp<_currentModel->numFloatingSpecies; isp++)
            {

                string speciesId = _currentModel->sp_list[isp].id;
                if(speciesId == innerString)
                {

                    string compartment = "vol__" + _currentModel->sp_list[isp].compartment;
                    bool isUnitVolume = _currentModel->compartmentsList[_currentModel->sp_list[isp].compartment] == 1.0;
                    if (divideVolumes)
                        replaceStream += "("

                    replaceStream += "x(" + str(isp+1) + ")"

                    if (divideVolumes):
                    {
                        if (!isUnitVolume):
                            replaceStream += "/" + compartment ;

                        replaceStream += ")";
                    }
                    match = True
                    break;
                }

            }

            if (!match )
            {
                if (innerString == "exponentiale")
                    replaceStream << "exp(1)";
                else if (innerString == "INF")
                    replaceStream << "Inf";
                else if (innerString == "arcsin")
                    replaceStream << "asin";
                else if (innerString == "arccos")
                    replaceStream << "acos";
                else if (innerString == "arctan")
                    replaceStream << "atan";
                else if (innerString == "arcsec")
                    replaceStream << "asec";
                else if (innerString == "arccsc")
                    replaceStream << "acsc";
                else if (innerString == "arccot")
                    replaceStream << "acot";
                else if (innerString == "arcsinh")
                    replaceStream << "asinh";
                else if (innerString == "arccosh")
                    replaceStream << "acosh";
                else if (innerString == "arctanh")
                    replaceStream << "atanh";
                else if (innerString == "arcsech")
                    replaceStream << "asech";
                else if (innerString == "arccsch")
                    replaceStream << "acsch";
                else if (innerString == "arccoth")
                    replaceStream << "acoth";

                else 
                    replaceStream << innerString;
            }
        }

        return replaceStream
    

def subConstants(equation, reactionId, speciesData, divideVolumes = True):
		
    # set up the scanner
   lx = Lexer(rules, skip_whitespace=True)
   lx.input('equation')

   try:
      for tok in lx.tokens():
          if tok.type == "IDENTIFIER":
             currentToken = tok.val						
             result << ReplaceStringToken(currentToken,  reactionId, divideVolumes);
    
          if tok.type == "FLOAT":
             result += tok.val 
                
          if tok.type == 'PLUS':
             result += '+'
            
          if tok.type == 'MINUS':
             result += '-'
            
          if tok.type == 'MULTIPLY':
             result += '*'
            
          if tok.type == 'MINUS':
             result += '-'
            
          if tok.type == 'DIVIDE':
             result += '/'
            
          if tok.type == 'POWER':
             result += '^'
            
          if tok.type == 'LP':
             result += '('
            
          if tok.type == 'RP':
             result += ')'
            
          if tok.type == 'COMMA':
             result += ','        
        
        result += ';'
        return result
    
   except LexerError as err:
      print('LexerError at position %s' % err.pos)


# subConstants function for matlab columns, for example, produces x(:,4) instead of x(4)
def subConstantsCol(const string &equation, const string &reactionId, bool divideVolumes = true)
		# set up the scanner
		stringstream ss(equation);

		TScanner scanner;
		scanner.setStream(&ss);
		scanner.startScanner();
		scanner.nextToken();

		stringstream result;

		try
		{
			while(scanner.getToken()!= tEndOfStreamToken)
			{

				switch ( scanner.getToken() )
				{
				case tWordToken :	
					{
						string currentToken = scanner.tokenToString ( scanner.getToken() );						
						result << ReplaceStringTokenCol(currentToken,  reactionId, divideVolumes);
					}
					break;
				case tDoubleToken:
					result << scanner.tokenDouble;
					break;
				case tIntToken	 :  
					result << scanner.tokenInteger;
					break;
				case tPlusToken	 :	result << "+";
					break;
				case tMinusToken :  result << "-";
					break;
				case tDivToken	 :  result << "/";
					break;
				case tMultToken	 :	result << "*";
					break;
				case tPowerToken :	result << "^";
					break;
				case tLParenToken:	result << "(";
					break;
				case tRParenToken:	result << ")";
					break;
				case tCommaToken :	result << ",";
					break;
				default			 :  MatlabError* ae =
										new MatlabError("Unknown token in subConstants (matlabTranslator): " + scanner.tokenToString( scanner.getToken() ));
					throw ae;
				}
				scanner.nextToken();
			}
			result << ";";
		}
		catch (MatlabError ae)
		{
			throw ae;
		}
		return result
	}

	// if expression is parenthesized, drop parentheses
	string stringInside(const string &str)
	{
		if (str == "" || str.length() < 2) return str;
		return str.substr(1, str.length()-2);
	}
