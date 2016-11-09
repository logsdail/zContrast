/*
 *  aMinimisationDetails.h
 *  zContrast
 *
 *  Created by Andrew Logsdail on 13/08/2010.
 *  Copyright 2010 University of Birmingham. All rights reserved.
 *
 */

#ifndef MINIMISATIONDETAILS_H
#define MINIMISATIONDETAILS_H
#include <vector>
#include <math.h>
#include <cstdlib>
#include <iostream>
#include "Structures.h"

class LocalMinimisationDetails{
	
public:
	LocalMinimisationDetails() 
	{ init(); }
	
	~LocalMinimisationDetails() 
	{;}
	
	void init();
	
	Direction getDirections();
	
	bool getBUni() 
	{ return bUni; }
	
	bool getUniFinished();
	
	void setBUni();
	
	bool getBMulti() 
	{ return bMulti; }
	
	void setBMulti();
	
	void setBPowell();
	
	bool getBPowell() 
	{ return bPowell; }
	
	//RotationPoint getLastPoint()
	//{ return lastPoint; }	
	
	//void setLastPoint(const RotationPoint &rp) 
	//{ lastPoint = rp; }
	
	// void nextDirection(const RotationPoint &rp);
	
	void setDirection();
	
	void setSearchLimits(const LinearStruct &ls) 
	{ searchLimits = ls; }
	
	void nextUnivariate(const int steps);

private:
	void setDefaults();
	
	void setSearchVectors(const float &dg);
	
	Direction minimisePowell();
	
	std::vector<Direction> search_vectors;
	
	std::vector<int> search_directions;
	
	bool bUni; // Univariate
	bool bMulti; // Multivariate	
	bool bPowell; //Powell
	
	//RotationPoint lastPoint;
	LinearStruct searchLimits;
	
	// Search Rotation //
	
	int total_vectors;
	int uni_current;
	int uni_count;
	int powell_degrees;
	
	////////////////////////
};
#endif
