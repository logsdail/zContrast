/*
 *  OptimisationUtils.cpp
 *  zContrast
 *
 *  Created by Andrew Logsdail on 02/06/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

/**
 02/06/2011
 - Removed parallelisation from random point selection
 **/

#include "GlobalOptimisationUtils.h"

using namespace std;

vector<RotationPoint> getRandomPoints(const int s, int *idum, LinearStruct ls)
// Sets RotationPoint variables randomly
// Returns vector of values;
{	
	RotationPoint rp;
	vector<RotationPoint> vRP;
	
	/**	
	 #pragma omp parallel for ordered default(none) \
	 private(rp) \
	 shared(vRP)
	 **/
	
	for (int i = 0; i < s; i++)
	{	
		rp = getRandomPoint(idum, ls);	
		// #pragma omp critical
		//		{
		vRP.push_back(rp);
		//		}
	}
	
	return vRP;
}

RotationPoint getRandomPoint(int *idum, LinearStruct ls)
// Sets RotationPoint variables randomly
// Inputs - size - number of points needed
//			ls - All variables associated with rotation maxima
{
	int randomNum;
	RotationPoint temp;
	
	randomNum = randomNumber((ls.theta_no_steps + 1),idum);	
	temp.theta = ls.theta_min + (randomNum * ls.theta_step);
	
	//randomNum = randomNumber((searchLimits.phi_no_steps + 1),idum);	
	//temp.phi = searchLimits.phi_min + (randomNum * searchLimits.phi_step);
	
	// This will distribute points evenly *hopefully* across the surface of a sphere.
	// Needs testing - should generate better convergence to minima
	
	temp.phi = ls.phi_max + 1;
	while ((temp.phi > ls.phi_max) || (temp.phi < ls.phi_min))
	{
		randomNum = (int) round(acos((randomNumberF(2,idum)) - 1) * RAD2DEG);
		randomNum = (int) round(randomNum / ls.phi_step);
		temp.phi = ls.phi_min + (randomNum * ls.phi_step);
	}
	
	randomNum = randomNumber((ls.psi_no_steps + 1),idum);	
	temp.psi = ls.psi_min + (randomNum * ls.psi_step);
	
	return tidy(temp);
}

int getMinimum(const vector<RotationPoint> &vec, const RotationPoint rp, vector<string> *o)
// Find minimum in vector
// Returns data position from that point in vector current
{	
	vector<RotationPoint> v = vec;
	int v_size = v.size();
	int low = 0;
	
	if (v_size == 0) 
	{
		emptyVectorError(o);
	}
	else
	{
		// If we have been given the value to use as the previous minimum, find this and compare to other values
		// Hopefully this shouldn't slow things down. I'll parallelise trivially to make sure.
		
		if (rp.value != -1)
		{
			float tx, ty, tz;
			
			// Parallelised for a little extra speed //
#pragma omp parallel for ordered default(none) \
shared(low,v_size,v) private(tx,ty,tz)
			for (int i = 0; i < v_size; i++)
			{
				tx = fabs(v[i].theta - rp.theta);
				ty = fabs(v[i].phi - rp.phi);
				tz = fabs(v[i].psi - rp.psi);
				
				
				if ((tx < ERRORBAR) && (ty < ERRORBAR) && (tz < ERRORBAR))
				{
					low = i;
				}
			}
		}
		
		// Parallelised for a little extra speed //
#pragma omp parallel for ordered default(none) \
shared(low,v_size,v)
		for (int i = 0; i < v_size; i++)
		{
			if (v[low].value > v[i].value)
			{
				low = i;
			}
		}
	}
	
	return low;
}

int getMaximum(const vector<RotationPoint> &vec, vector<string> *o)
// Find minimum in vector
// Returns data position from that point in vector current
{
	vector<RotationPoint> v = vec;
	int v_size = v.size();
	int high = 0;
	
	if (v_size == 0) 
	{
		emptyVectorError(o);
	}
	else
	{
		// Parallelised for a little extra speed //
#pragma omp parallel for ordered default(none) \
shared(high,v_size,v)
		for (int i = high + 1; i < v_size; i++)
		{
			if (v[high].value < v[i].value) 
			{
				high = i;
			}
		}
	}		
	
	return high;
}

RotationPoint emptyRotationPoint()
// Return empty RP
{ 
	RotationPoint r;
	r.theta = -1;
	r.phi = -1;
	r.psi = -1;
	r.value = -1; 
	return r; 
}

RotationPoint tidy(RotationPoint rp)
// Get rid of rounding errors from our value.
// Only really adds cosmetic value
// Inputs: RotationPoint - To be tidied
// Returns: RotationPoint - Tidied
{
	
	// Check modulus of remainder //
	float error = fabs(fmod(rp.theta,1));
	
	if (error < ERRORBAR || (1-error) < ERRORBAR)
		rp.theta = round(rp.theta);
	
	//Y////////////////////
	
	// Check modulus of remainder //
	error = fabs(fmod(rp.phi,1));
	
	if (error < ERRORBAR || (1-error) < ERRORBAR)
		rp.phi = round(rp.phi);
	
	//Z////////////////////
	
	// Check modulus of remainder //
	error = fabs(fmod(rp.psi,1));
	
	if (error < ERRORBAR || (1-error) < ERRORBAR)
		rp.psi = round(rp.psi);
	
	///////////////////////
	
	return rp;
}

void emptyVectorError(vector<string> *o)
// Displays error that we are going to check an empty vector for values
{
	o->push_back("!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!");
	o->push_back("!! Vector has no entries so we will be          !!");
	o->push_back("!! seeing a segmentation fault. Did you specify !!");
	o->push_back("!! a comparative value?                         !!");
	//o->push_back("!! Currently : " + getReference());
	o->push_back("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
}

