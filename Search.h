#ifndef SEARCH_H
#define SEARCH_H

/**
 02/06/2011
		- Commented out headers inherited from Utils.h
 **/

// #include <cstdlib>
// #include <iostream>
#include "CVariables.h"
// #include "Structures.h"
#include "CreateRotations.h"
// #include "Utils.h"
// #include "GlobalOptimisationUtils.h"
#include "GlobalOptimisation.h"
#include "LocalMinimisationDetails.h"

class Search{
	
public:
	Search() {output_content = NULL;}
	
	Search(int *s, std::vector<std::string> *o)
	{ init(s,o); }
	
	~Search() 
	{ clearUp(); }
	
	void init(int *s, std::vector<std::string> *o);
	
	void clearUp();
	
	void getVariables(SearchV seV, const Variables &v);
	
	bool getOptimisation() 
	{ return bOptimisation; }
	
	bool getGlobal() 
	{ return bGlobal; }
	
	bool getGlobalConverged() 
	{ return optimisation.getConverged(); }
	
	std::string getOptimisationRef() 
	{ return optimisation.getReference(); }
	
	// Set value for optimisation point
	// Inputs (x) - X coordinate
	//        (y) - Y coordinate
	//        (z) - Z coordinate
	//	  (value) - Value
	void setOptimisationValue(const float &x, const float &y, const float &z, const float &value)
	{ RotationPoint temp; temp.theta = x; temp.phi = y; temp.psi = z; temp.value = value; current_results.push_back(temp); }
	
	bool getComplete();
	
	Variables getNext(const unsigned int current);
	
	int getRotationsArraySize() 
	{ return rotations.size(); }
	
	void setComplete() 
	{ bComplete = true; }
	
	void setIncomplete() 
	{ bComplete = false; }
	
	void setVariables(const Variables v) 
	{ var = v; } 
	
	// RotationPoint getCurrentMin()
	// { return optimisation.getCurrentMin(); }
	
	int getHistoryCopyCounter()
	{ return historyCopyCounter; }
	
	RotationPoint getMinimumPoint()
	{return optimisation.getMinimumPoint();}

private:
	bool bLinear;
	bool bGlobal;
	bool bOptimisation;
	bool bScreen;
	bool bPeriodic;
	bool bOptimised;
	bool bConverged;
	bool bComplete;
	bool bOptimising;
	bool bInTheLoop; // Using this to identify if we are looping
	
	std::vector<RotationPoint> rotations;
	std::vector<RotationPoint> possible_points;
	std::vector<RotationPoint> current_results;
	std::vector< std::vector< std::vector<float> > > historyNew; // 3D Matrix of results
	std::vector<std::string> *output_content;
	
	GlobalOptimisation optimisation;
	Variables var;
	LocalMinimisationDetails miniDetails;
	
	int historyCopyCounter;
	int powellStepCounter;
	
	void setStartingPoints(const LinearStruct linearStruct);
	
	void setOptimisationStartingPoints();
	
	void addCurrentToHistory();

	void copyFromHistory();
	
	void locallyMinimised(const int &i);
	
	void printHistoryArrayError(RotationPoint rp);
};
#endif
