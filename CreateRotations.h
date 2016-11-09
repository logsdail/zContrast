#ifndef CREATE_ROTATIONS_H
#define CREATE_ROTATIONS_H

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 23/08/2011
 - Edited pbc checks to incorporate univariate searches
 - Created pbc_check method
 **/

// #include <cstdlib>
// #include <iostream>
// #include <vector>
// #include <math.h>
// #include <omp.h>
// #include "Structures.h"
// #include "Utils.h"
#include "GlobalOptimisationUtils.h"

class Create_Rotations{
	
public:
	Create_Rotations() {;}
	~Create_Rotations() {;}
	
	std::vector<RotationPoint> createPoints(const bool &bMulti, const LinearStruct &linearStruct, const LinearStruct &limits, const bool &pbc); 	

	std::vector<RotationPoint> createPoints(const LinearStruct &linearStruct);
	
private:
	std::vector<RotationPoint> createPointsLinear(const LinearStruct linearStruct);
	
	std::vector<RotationPoint> createPointsUnivariate(const LinearStruct &linearStruct);
	
	std::vector<RotationPoint> pbc_check(const bool pbc, std::vector<RotationPoint> rotations, const LinearStruct linearStruct, const LinearStruct &limits);
};
#endif

