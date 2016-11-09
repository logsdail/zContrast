/*
 *  readVariables.h
 *  zContrast
 *
 *  Created by Andrew Logsdail on 20/01/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

#ifndef READVARIABLES_H
#define READVARIABLES_H

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 **/

// #include <cstdlib>
#include <iostream>
#include <fstream>
// #include "Structures.h"
#include "Utils.h"

class ReadVariables
{
public:
	
	ReadVariables();
	
	~ReadVariables() 
	{;}
	
	bool openFile(const std::string &filename);
	
	CreateV getCreateVariables() 
	{return cr_v;}
	
	CompareV getCompareVariables() 
	{return co_v;}
	
	SearchV getSearchVariables() 
	{return se_v;}
	
	GaV getGaVariables() 
	{return se_v.ga_v;}

	// Variable carriers
	MainV ma_v;
	CreateV cr_v;
	CompareV co_v;
	SearchV se_v;
	
private:
	bool getVariables(std::ifstream &inData);
	
	bool mainVariable(const std::string variable, std::string value);
	
	bool createVariable(const std::string variable, std::string value);
	
	bool compareVariable(const std::string variable, std::string value);
	
	bool searchVariable(const std::string variable, std::string value);
	
	SearchV checkOrder(SearchV seV);
	
	bool gaVariable(const std::string variable, std::string value);
	
	bool readParameters();
	
};
#endif
