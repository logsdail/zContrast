#ifndef OPTIMISATIONMUTATION_H
#define OPTIMISATIONMUTATION_H

/*
 *  OptimisationMutation.h
 *  zContrast
 *
 *  Created by Andrew Logsdail on 02/06/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 **/

// #include <cstdlib>
// #include <iostream>
// #include <vector>
// #include <string>
// #include "Structures.h"
// #include "Utils.h"
#include "GlobalOptimisationUtils.h"

class OptimisationMutation{
	
public:
	
	OptimisationMutation() {idum = NULL;}
	~OptimisationMutation() {;}
	
	void setVariables(int *s, std::vector<std::string> *o, bool b, const std::string mutation);
	
	std::vector<RotationPoint> getMutantPoints(int numberOfMutants, int current_generation_count,
											   LinearStruct ls,
											   const std::vector<RotationPoint> points_copy);
	
private:
	
	float newRotationAngle(const float current_value, const int no_steps, const float step_size, 
						   const float min, const float max, const float range, const int dynamic_breadth);
	
	int *idum; // Pointer for random number generator
	std::vector<std::string> *output_content;
	
	// MUTATION TYPES //
	bool bMStaticPoint;
	bool bMDynamicPoint;
	bool bScreen;
};
#endif

