#ifndef OPTIMISATION_H
#define OPTIMISATION_H

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 **/

// #include <cstdlib>
// #include <iostream>
// #include <vector>
// #include <string>
// #include "Structures.h"
// #include <omp.h>
// #include "math.h"
// #include "Utils.h"
#include "GlobalOptimisationMutation.h"
#include "GlobalOptimisationOffspring.h"
// #include "GlobalOptimisationUtils.h"

class GlobalOptimisation{
	
public:
	GlobalOptimisation() {idum = NULL; output_content = NULL;}
	
	~GlobalOptimisation() {;}
	
	GlobalOptimisation(int *s, std::vector<std::string> *o) 
	{ init(s,o); }
	
	void init(int *s, std::vector<std::string> *o);

	void setTempPoint(const RotationPoint &l) 
	{ tempPoint = l; }
	
	void setTempPoint(const float t, const float ph, const float ps, const float v)
	{ tempPoint.theta = t; tempPoint.phi = ph; tempPoint.psi = ps; tempPoint.value = v;}
	
	void setPoint(const RotationPoint &l);
	
	void setPoints(const std::vector<RotationPoint> &rp);
	
	// Sets random start points
	void setRandomStartPoints() 
	{newPoints = getRandomPoints(size,idum,getSearchLimits());}
	
	void setReference(const std::string &r)
	{ reference = r; }
	
	void setSearchLimits(const LinearStruct &sl) 
	{ searchLimits = sl; }
	
	void setSize(const int &s) 
	{ size = s; points.resize(s); }

	RotationPoint getTempPoint() 
	{ return tempPoint; }
	
	RotationPoint getNewPoint();
	
	std::vector<RotationPoint> getNewPoints();
	
	LinearStruct getSearchLimits() 
	{ return searchLimits; }
	
	int getSize() 
	{ return size; }
	
	bool getPopulationFull();
	
	int getPopulationSize() 
	{ return size; }
	
	bool getMutantsFull();
	int getMutantsSize() 
	{ return numberOfMutants; }
	
	bool getOffspringFull();
	int getOffspringSize() 
	{ return numberOfOffspring; }
	
	bool getConverged() 
	{ return bConverged; }
	
	bool getbScreen() 
	{ return bScreen; }
	
	std::string getReference()
	{ return reference; }

	void newPopulation();	
	
	RotationPoint getCurrentMin()
	{ return minimum; }
	
	LinearStruct minimisationPoints(const Direction &d);
	
	void setVariables(GaV ga_v);
	
	void printMinimum();
	
	void printTempPoint(const std::string &front);
	
	RotationPoint getMinimumPoint()
	{return minimum;}
	
	void setMinimumPoint(RotationPoint rp)
	{minimum = rp;}
	
private:
	void setUp();
	
	void resetPoints();
	
	void checkConverged();
	
	void next() 
	{ current++; }
	
	void nextMutant() 
	{ current_mutant++; }
	
	void nextOffspring() 
	{ current_offspring++; }
	
	void printSetPoints();
	
	std::string reference;
	
	OptimisationOffspring oo;
	
	OptimisationMutation om;
	
	std::vector<RotationPoint> points;
	std::vector<RotationPoint> newPoints;
	std::vector<RotationPoint> mutants;
	std::vector<RotationPoint> offspring;
	
	RotationPoint tempPoint;
	RotationPoint minimum;
	
	LinearStruct searchLimits;
	
	int size; // NUMBER IN OPTIMIZATION
	int current; // CURRENT NUMBER
	int generations; // NUMBER OF GENERATIONS BEFORE TERMINATION
	int current_generation_count; // COUNT OFF CURRENT ROTATIONS AT THIS LOW
	int tsize; // TOURNAMENT SIZE

	float mrate; // MUTATION RATE
	float noff; // PERCENTAGE OFFSPRING
	
	int numberOfMutants;
	int current_mutant;
	int numberOfOffspring;
	int current_offspring;	
	int *idum; // Pointer for random number generator
	std::vector<std::string> *output_content;
	
	bool bPopulationFull;
	bool bMutantsFull;
	bool bOffspringFull;
	bool bConverged;
	bool bScreen;
	// TYPES OF MUTATION
	bool bMChildren;
};
#endif
