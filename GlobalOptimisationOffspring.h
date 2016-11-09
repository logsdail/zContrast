#ifndef OPTIMISATIONOFFSPRING_H
#define OPTIMISATIONOFFSPRING_H

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
 - Imported Fitness function into this file for Roulette method
 **/

// #include <cstdlib>
// #include <iostream>
// #include <vector>
// #include <string>
// #include "Structures.h"
// #include "Utils.h"
#include "GlobalOptimisationUtils.h"

class OptimisationOffspring{
	
public:
	
	OptimisationOffspring() { tsize = 2; idum = NULL;}
	~OptimisationOffspring() {;}
	
	void setVariables(int *s, std::vector<std::string> *o, bool b, int tournament_size, 
					  const std::string type, const std::string parents, const std::string fitness);
	
	std::vector<RotationPoint> getOffspringPoints(int numberOfOffspring,
												  LinearStruct limits,
												  const std::vector<RotationPoint> points_copy);
	
	float fitness(const int rp, const std::vector<RotationPoint> points);
	
private:
	
	int *idum; // Pointer for random number generator
	std::vector<std::string> *output_content; // Outputs
	
	// MATING TYPES //
	bool bOUniformCrossover;
	bool bOTournament;
	bool bORoulette;
	bool bScreen;
	
	// FITNESS TESTS //
	bool bFExponential;
	bool bFLinear;
	bool bFTanh;
	
	int tsize; //Tournament
};
#endif

