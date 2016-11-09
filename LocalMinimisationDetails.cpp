/*
 *  aMinimisationDetails.cpp
 *  zContrast
 *
 *  Created by Andrew Logsdail on 13/08/2010.
 *  Copyright 2010 University of Birmingham. All rights reserved.
 *
 */

/**
 Updates:
 20/09/2011
 - Removed calculation of number of steps for Powell's
 - This is passed in now to save calculation.
 21/09/2011
 - Added method to setTempPoint for local minimisation routines
 **/

#include "LocalMinimisationDetails.h"

using namespace std;

void LocalMinimisationDetails::init()
// Set default values
{
	///// INITIALISATIONS /////
	//lastPoint.theta = 88888;
	//lastPoint.phi = 88888;
	//lastPoint.psi = 88888;
	//lastPoint.value = 88888;
	
	total_vectors = 0;
	uni_count = 0;
	uni_current = 0;
	///////////////////////////
}

void LocalMinimisationDetails::setBMulti()
// Set up for multivariate calculation
{
	bMulti = true;
	bUni = false;
	bPowell = false;
	setDirection();
}

void LocalMinimisationDetails::setBUni()
// Set up for univariate calculation
{
	bMulti = false;
	bPowell = false;
	bUni = true;
	setDirection();
}

void LocalMinimisationDetails::setBPowell()
// Set up for minimisation via Powell's method
{
	bMulti = false;
	bUni = true;
	bPowell = true;
	setDirection();
}

void LocalMinimisationDetails::setDirection() 
// Sets degrees between each univariate search
// Updated to just be three perpendicular axes x, y and z
{ 	
	total_vectors = 3; // x, y and z
	// total_vectors++;
	
	search_vectors.resize(total_vectors);
	search_directions.resize(total_vectors);
	
	if (bMulti)
	{
		for (int i = 0; i < total_vectors; i++)
		{
			search_vectors[i].theta = search_vectors[i].phi = search_vectors[i].psi = 1;
		}
	}
	else
	{
		// X AXES //
		search_vectors[0].theta = 1;
		search_vectors[0].phi = 0;
		search_vectors[0].psi = 0;
		// Y AXES //
		search_vectors[1].phi = 1;
		search_vectors[1].theta = 0;
		search_vectors[1].psi = 0;
		// Z AXES //
		search_vectors[2].psi = 1;
		search_vectors[2].phi = 0;
		search_vectors[2].theta = 0;
	}
	
	if (bPowell)
	{
		total_vectors++;
	}
}

bool LocalMinimisationDetails::getUniFinished()
// Check if search is finished for univariate minimisation
// Returns boolean yes or no
{
	if (uni_count == total_vectors)
	{
		//cout << "Yes: " << uni_count << endl;
		uni_count = uni_current = 0;
		return true;
	}
	
	//cout << "No: " << uni_count << endl;
	return false;
}

void LocalMinimisationDetails::nextUnivariate(const int steps)
// Gives new direction for univariate calculation
// Inputs(rp) - Last rotation point minima
{	
	//cout << diffx << " " << diffy << " " << diffz;
	//cout << steps << endl;
	
	// Create Total Displacement for Powell's Search Vector //
	if (bPowell && (uni_current < (total_vectors-1)))
	{
		search_directions[uni_current] = steps;
	}
	
	// CHECK IF WE ARE AT THE SAME POINT AS LAST TIME
	if (steps == 0)
	{
		uni_count++;
		//setLastPoint(rp);
	}
	else
	{
		uni_count = 0;
	}
	
	uni_current++;
	if (uni_current == total_vectors)
	{
		uni_current = 0;
	}
}
 
Direction LocalMinimisationDetails::getDirections()
// Method to get directions vectors
// Returns Direction variable
{
	if ((!bPowell) || ((bPowell) && (uni_current < (total_vectors-1))))
	{
		return search_vectors[uni_current];
	}
	else
	{
		return minimisePowell();
	}
}

Direction LocalMinimisationDetails::minimisePowell()
// Calculates cumulative search vector, and minimises it if possible
// Returns Direction Variable
{
	int i = 0;
	int max = 0;
	
	// FIND OUT MAXIMUM SEARCH VECTOR //
	for (i = 1; i < (total_vectors - 1); i++)
	{
		if (fabs(search_directions[i]) > search_directions[max])
		{
			max = i;
		}
	}

	int remain;
	int j;
	
    //cout << "Search Directions:" << search_directions[0] << " " << search_directions[1] << " " << search_directions[2] << endl;
	
	// DIVIDE BY ALL SEARCH VECTORS TO SEE IF FITS to BE MINIMISED//
	for (i = 2; i <= fabs(search_directions[max]); i++)
	{
		remain = 0;
		for (j = 0; j < (total_vectors - 1); j++)
		{
			remain += search_directions[j]%i;
		}
		
		// IF REMAINDER IS ZERO WE HAVE SUCCESS //
		if (remain == 0)
		{
			for (j = 0; j < (total_vectors - 1); j++)
			{
				search_directions[j] /= i;
			}
		}
	}
	
	// New direction
	Direction d;
	d.theta = 0;
	d.phi = 0;
	d.psi = 0;
	////////////////
	
	// cout << "Search Directions:" << search_directions[0] << " " << search_directions[1] << " " << search_directions[2] << endl;
	
	for (int i = 0; i < (total_vectors - 1); i++)
	{
		d.theta += search_directions[i]*search_vectors[i].theta;
		d.phi += search_directions[i]*search_vectors[i].phi;
		d.psi += search_directions[i]*search_vectors[i].psi;
	}
	
	//cout << "Powell's Vector:" << d.theta << " " << d.phi << " " << d.psi << endl;

	return d;
}
