/*
 *  OptimisationMutation.cpp
 *  zContrast
 *
 *  Created by Andrew Logsdail on 02/06/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

/**
 02/06/2011 - Outsourced mutation from Optimisation.cpp
 - Needs breaking down into smaller classes
 - So we have separated out newRotationAngles. Bit more readable now!
 **/

#include "GlobalOptimisationMutation.h"

using namespace std;

void OptimisationMutation::setVariables(int *s, vector<string> *o, bool b, const string mutation)
/**
 Method to set variables for mutation.
 Inputs: Seed for Random Number
 String - static or dynamic
 Print to Screen
 **/
{
	// Initiate Random Numbers
	idum = s;	
	output_content = o;
	//////////////////////////
	bMStaticPoint = false;
	bMDynamicPoint = false;
	
	// Assign values now
	bScreen = b;
	// MUTATION TYPES //
	if (cmpStr(mutation,"static"))
	{
		bMStaticPoint = true;
	}
	else if (cmpStr(mutation,"dynamic")) 
	{
		bMDynamicPoint = true;
	}
	else
	{
#pragma omp critical
		{
			output_content->push_back("Invalid mutation type: " + mutation);
		}
		//return EXIT_FAILURE;
	}
	
	//return EXIT_SUCCESS;
}

vector<RotationPoint> OptimisationMutation::getMutantPoints(int numberOfMutants, int current_generation_count,
															LinearStruct ls,
															const vector<RotationPoint> points_copy)
// Sets RotationPoint variables to mutated values
// Inputs:	- Number of Mutants needed
//			- Current generation count for dynamic mutation
//			- Linearstruct with possible ranges of values
//			- Copy of points to mutate
{
	if (bScreen) 
	{
		string mutantS;
		string numMutantsS;
		int mutantI = points_copy.size();
		NumberToString(mutantI,mutantS);
		NumberToString(numberOfMutants,numMutantsS);
		
#pragma omp critical
		{
			output_content->push_back("Working out mutants : " + numMutantsS + ". Points size is : " + mutantS);
		}
	}
	
	// New mutants array. We'll return this at the end
	vector<RotationPoint> new_mutants(numberOfMutants);
	
	if ((!bMStaticPoint) && (!bMDynamicPoint))
	{
#pragma omp critical
		{
			output_content->push_back("Error: No Mutation Scheme defined");
		}
		new_mutants = getRandomPoints(numberOfMutants,idum,ls);
	}
	else
	{
		// Set to one, as this is the value for static mutation
		int breadth = 1;
		
		if (bMDynamicPoint) 
		{
			breadth = current_generation_count;
		}
		/**
		 // Have had trouble with this parallelisation //
		 #pragma omp parallel for ordered default(none) \
		 shared(ls, cout, breadth, points_copy, new_mutants)
		 **/
		for (int i = 0; i < numberOfMutants; i++)
		{
			int option = randomNumber((int) 3,idum); // Get random number to represent one of the three variables
			int position = randomNumber(points_copy.size(),idum);
			RotationPoint rp = points_copy[position];
			
			// This has been edited as the arrays are passed in early now
			// Instead of editing mutants directly.
			// position = randomNumber(points_copy.size(),idum); // Get random point
			
			/**			#pragma omp critical
			 {	**/
			//rp = points_copy[position];
			/**			}	**/
			
			if (option == 0)
				// Change x variable
			{
				rp.theta = newRotationAngle(rp.theta,ls.theta_no_steps,ls.theta_step, ls.theta_min, ls.theta_max, ls.theta_range, breadth);
			}
			else if (option == 1)
				// Change y variable
			{
				rp.phi = newRotationAngle(rp.phi,ls.phi_no_steps,ls.phi_step, ls.phi_min, ls.phi_max, ls.phi_range, breadth);
			}
			else if (option == 2)
				// Change z variable
			{
				rp.psi = newRotationAngle(rp.psi,ls.psi_no_steps,ls.psi_step, ls.psi_min, ls.psi_max, ls.psi_range, breadth);
			}
			
			// Set current value to zer0 to prevent any odd values appearing
			rp.value = 0;
			new_mutants[i] = rp;
		}
	}
	
	if (bScreen) 
	{
		string mutantS;
		int mutantI = new_mutants.size();
		NumberToString(mutantI,mutantS);
		
#pragma omp critical
		{
			output_content->push_back("Mutants generated : " + mutantS);
		}
	}
	return new_mutants;
}

float OptimisationMutation::newRotationAngle(const float current_value, const int no_steps, const float step_size, 
											 const float min, const float max, const float range, const int dynamic_breadth)
/**
 Calculates new possible values for mutants
 Inputs:	- current_value is the current rotation angle
 - number of possible steps we can choose from
 - step size of these steps
 - minimum value
 - maximum value
 - minimum - maximum i.e. range
 - dynamic_breadth
 **/
{
	vector<float> newValues; // Sometimes appears empty? Comparison problem.
	
	// Changed back to cater for square grid
	for (signed int tx = 0; tx <= no_steps; tx++ )
	{
		newValues.push_back(min + (tx * step_size));
	}
	// Now we have all possible new values
	
	// CREATE DYNAMIC SOLUTION //
	// Here the possible mutation vectors decrease in size with the generations we have been at this minima //
	if (dynamic_breadth > 1)
	{
		float spread_max = ((range + 1) / (2 * dynamic_breadth));
		
		for (unsigned int j = 0; j < newValues.size(); j++)
			// Check if we are beyond limits
			if ((newValues[j] > (current_value + spread_max)) || (newValues[j] < (current_value - spread_max)))
			{
				newValues.erase(newValues.begin()+j);
				j--;
			}
	}
	///////////////////////////////
	
	if (newValues.size() > 0) 
	{
		return newValues[randomNumber(newValues.size(),idum)];
	}
	else 
	{
		string current_valueS;
		string dynamic_breadthS;
		NumberToString(current_value,current_valueS);
		NumberToString(dynamic_breadth,dynamic_breadthS);
		
#pragma omp critical
		{
			output_content->push_back("");
			output_content->push_back("Minor error : Mutation Unsuccessful for " + current_valueS + ". Current breadth = " + dynamic_breadthS);
			output_content->push_back("");
		}
		
		return current_value;
	}
	
}
