/*
 *  OptimisationUtils.h
 *  zContrast
 *
 *  Created by Andrew Logsdail on 02/06/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

#ifndef __OPTIMISATIONUTILS_HPP__
#define __OPTIMISATIONUTILS_HPP__

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 23/08/2011
 - Added in tidy function for Rotation Points
 **/

// #include <cstdlib>
#include <iostream>
// #include <vector>
#include "omp.h"
#include "Utils.h"
// #include "Structures.h"

// Random Points

std::vector<RotationPoint> getRandomPoints(const int s, int *idum, LinearStruct ls);

RotationPoint getRandomPoint(int *idum, LinearStruct ls);

// Minima and Maxima of arrays
int getMinimum(const std::vector<RotationPoint> &vec, const RotationPoint rp, std::vector<std::string> *o);

int getMaximum(const std::vector<RotationPoint> &vec, std::vector<std::string> *o);

// Tidy up Rotation Point
RotationPoint tidy(RotationPoint rp);

// Initialise empty rotation
RotationPoint emptyRotationPoint();

void emptyVectorError(std::vector<std::string> *o);


#endif

