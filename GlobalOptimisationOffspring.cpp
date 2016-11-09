/*
 *  OptimisationOffspring.cpp
 *  zContrast
 *
 *  Created by Andrew Logsdail on 02/06/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

/**
 02/06/2011 - Outsourced offspring from Optimisation.cpp
 - Needs breaking down into smaller classes
 - We have imported fitness calculation into this class
 **/

#include "GlobalOptimisationOffspring.h"

using namespace std;

void OptimisationOffspring::setVariables(int *s, vector<string> *o, bool b, int tournament_size, 
										 const string type, const string parents, const string fitness)
/**
 Method to set variables for mutation.
 Inputs: Seed for Random Number
 String - Mating Type
 String - Parent Selection
 String - Fitness Type
 Print to Screen
 **/
{
	// Initiate Random Numbers
	idum = s;	
	output_content = o;
	// Set Tournament size
	tsize = tournament_size;
	//////////////////////////
	bOUniformCrossover = false;
	bOTournament = false;
	bORoulette = false;
	bFExponential = false;
	bFLinear = false;
	bFTanh = false;
	
	// Assign values now
	bScreen = b;
	// CROSSOVER //
	if (cmpStr(type,"uniform")) 
	{
		bOUniformCrossover = true;
	}
	else
	{
#pragma omp critical
		{
			output_content->push_back("Invalid mating type: " + type);
		}
		//return EXIT_FAILURE;
	}
	
	// CROSSOVER SELECTION //
	if (cmpStr(parents,"tournament")) 
	{
		bOTournament = true;
	}
	else if (cmpStr(parents,"roulette")) 
	{
		bORoulette = true;
	}
	else
	{
#pragma omp critical
		{
			output_content->push_back("Invalid crossover selection type: " + parents);
		}
		//return EXIT_FAILURE;
	}	
	
	// FITNESS TESTS //
	if (cmpStr(fitness,"exponential")) 
	{
		bFExponential = true;
	}
	else if (cmpStr(fitness,"linear"))
	{
		bFLinear = true;
	}
	else if (cmpStr(fitness,"tanh"))
	{
		bFTanh = true;
	}
	else
	{
#pragma omp critical
		{
			output_content->push_back("Invalid fitness type: " + fitness);
		}
		//return EXIT_FAILURE;
	}
	
	//return EXIT_SUCCESS;
}

vector<RotationPoint> OptimisationOffspring::getOffspringPoints(int numberOfOffspring,
																LinearStruct limits,
																const std::vector<RotationPoint> points_copy)
// Sets RotationPoint variables to offspring values
{
	if (bScreen)
	{
		string snum;
		string spoints;
		int ipoints = points_copy.size();
		NumberToString(numberOfOffspring,snum);
		NumberToString(ipoints,spoints);
		
#pragma omp critical
		{
			output_content->push_back("Working out offspring : " + snum + ". Points size is : " + spoints);
		}
	}
	// Looks like the problem we are having is "points" is being overwritten or emptied mid run?
	// So what we'll do is create a copy of it just for this method - now takes points in
	vector<RotationPoint> new_offspring;
	// Let's see if this clears up problems. We'll return a new copy off the offspring to be assigned
	
	// offspring.clear();
	if (!bOUniformCrossover) 
	{
#pragma omp critical
		{
			output_content->push_back("Error: No mating scheme defined");
		}
		new_offspring = getRandomPoints(numberOfOffspring,idum,limits);
	}
	else 
	{
		int one = 0;
		int two= 0; 
		int option = 0;
		RotationPoint rp = emptyRotationPoint();
		/**	
		 // Parallel loop to generate uniform crossover //
		 #pragma omp parallel for default(none) \
		 private(one,two,option,rp) shared(points_copy,new_offspring) 
		 **/		
		for (int i = 0; i < numberOfOffspring; i++)
		{
			one = two = 1;
			// Make sure we don't get the same parents //
			while (one == two)
			{
				// Hopefully this will deal with infinite loops
				// Going to remove bottom else loop as a result
				one = randomNumber(points_copy.size(),idum);
				two = randomNumber(points_copy.size(),idum);
				
				// Select parents by tournament method // 
				if (bOTournament)
				{
					vector<RotationPoint> tourn;
					vector<int> positions;
					
					for (int j = 0; j < tsize; j++)
					{
						int pos = randomNumber(points_copy.size(),idum);
						positions.push_back(pos);
						
						/**						#pragma omp critical
						 { **/
						tourn.push_back(points_copy[pos]);
						/**	} **/
					}
					
					int t = getMinimum(tourn,emptyRotationPoint(),output_content);
					one = positions[t];
					positions.erase(positions.begin()+t);
					tourn.erase(tourn.begin()+t);
					
					t = getMinimum(tourn,emptyRotationPoint(),output_content);
					two = positions[t];
					positions.clear();
					tourn.clear();
				}
				else if (bORoulette)
					// #pragma omp critical
					// Select parents via roulette method instead of tournament //
				{
					// int counter = 0;				// We are going to use this as an escape route
					// int counter_limit = 1000;		// For this weird infinite loop we are seeing
					
					float compareF = randomNumber(idum);
					string vS;
					string fitnessS;
					string compareS;
					float fitnessF = fitness(one,points_copy);
					
					while (fitnessF >= compareF)
					{ 
						NumberToString(one,vS);
						NumberToString(fitnessF,fitnessS);
						NumberToString(compareF,compareS);
						
#pragma omp critical
						{
							output_content->push_back("Roulette failed (1) : " + vS + " : " + fitnessS + " > " + compareS);
						}
						
						one = randomNumber(points_copy.size(),idum);
						fitnessF = fitness(one,points_copy);
						compareF = randomNumber(idum);
						
						/** We've had problems with this infinitely looping **/
						// counter++;
						// if (counter > counter_limit)
						// {
						//	cout << "Exceeded roulette limit (1), breaking out. Points size is " << points_copy.size() << endl;
						//	compareF = 1;
						// }
						/** So this gives us a break out **/
					}
					
					// Reset counter
					// counter = 0;
					
					fitnessF = fitness(two,points_copy);
					compareF = randomNumber(idum);
					while (fitnessF >= compareF) 
					{
						NumberToString(two,vS);
						NumberToString(fitnessF,fitnessS);
						NumberToString(compareF,compareS);
#pragma omp critical
						{
							output_content->push_back("Roulette failed (2) : " + vS + " : " + fitnessS + " > " + compareS);
						}
						
						two = randomNumber(points_copy.size(),idum);
						fitnessF = fitness(two,points_copy);
						compareF = randomNumber(idum);
						
						/** We've had problems with this infinitely looping **/
						// counter++;
						// if (counter > counter_limit)
						// {
						// 	cout << "Exceeded roulette limit (2), breaking out. Points size is " << points_copy.size() << endl;
						//	compareF = 1;
						// }
						/** So this gives us a break out **/
					}
				}
			}
			
			// Create children if the parents aren't the same //
			
			/**			#pragma omp critical
			 {	**/
			rp = points_copy[one];
			
			int option_total = 0;
			
			for (int j = 0; j < 3; j++)
			{
				option = randomNumber((int) 2,idum);
				
				if ((j == 0) && (option == 1)) 
				{
					rp.theta = points_copy[two].theta;
				}
				else if ((j == 1) && (option == 1)) 
				{
					rp.phi = points_copy[two].phi;
				}
				else if ((j == 2) && (option == 1)) 
				{
					rp.psi = points_copy[two].psi;
				}
				
				option_total += option;
			}
			
			// Set value. If previously calculated set as so.
			if (option_total == 3) 
			{
				rp.value = points_copy[two].value;
			}
			else if (option_total != 0) 
			{
				rp.value = 88888;
			}
			// Push back new point //
			/**			}	**/
			
			new_offspring.push_back(rp);
		}
	}
	
	if (bScreen) 
	{
		string offspringS;
		int offspringI = new_offspring.size();
		NumberToString(offspringI,offspringS);
		
#pragma omp critical
		{
			output_content->push_back("Offspring generated : " + offspringS);
		}
	}
	
	return new_offspring;
}

float OptimisationOffspring::fitness(const int rp, const vector<RotationPoint> points)
// Work out fitness of a value compared to overall population
// Input: Current Rotation Point Index
//		: Current Selection of Points
// Output : float of fitness value
{
	int worst = getMaximum(points,output_content); // MAX
	int best = getMinimum(points,emptyRotationPoint(),output_content); // MIN
	float p = ((points[rp].value - points[best].value)/(points[worst].value - points[best].value)); // Normalise
	// EXPONENTIAL FACTOR //
	const float alpha = 3;
	// Should we soft code this? Probably //
	////////////////////////
	
	if (bFLinear) 
	{
		return (1-(0.7*p)); // Linear function
	}
	else if (bFExponential) 
	{
		return exp(-alpha*p); // Exponential function
	}
	else if (bFTanh) 
	{
		return 0.5*(1-tanh((2*p)-1)); // Hyperbolic Tangent function
	}
	else 
	{
		return p;
	}
}
