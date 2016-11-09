#include "GlobalOptimisation.h"

/** Updates:
 30/05/2011 
 - Added couts to offspring methods. Both Tournament and Roulette sometimes give infinite loop?
 - Changed setOffspringPoints and setMutantPoints to take the arrays as inputs, as opposed to
 editing the arrays mutants and offspring directly
 31/05/2011	
 - Remove parallelisation from setOffspringPoints and setMutantPoints. Continue to get segmentation
 faults from accessing arrays. This could be solved in the long run by outsourcing these files
 to separate classes and working on the access techniques
 02/06/2011
 - Outsource mutation methods to GlobalOptimisationMutation
 - Outsourced random points method to GlobalOptimisationUtils
 - Outsourced offsrping methods to GlobalOptimisationOffspring, along with fitness tests
 - Changed GA population selection to just elitist.
 **/


using namespace std;

void GlobalOptimisation::init(int *s, vector<string> *o) 
// Constructor
// Input - Seed for random number
{
	// Initiate Random Numbers
	idum = s;	
	output_content = o;
	//////////////////////////
	setUp();
}

void GlobalOptimisation::setUp()
// Setup class for fresh run
// Takes no inputs but configures for no mess
{
	bConverged = false;
	
	// GA PARAMETERS //
	// TO SCREEN //
	bScreen = false;	
	// MUTATION TYPES //
	bMChildren = false;
	
	minimum.theta = LARGE;
	minimum.phi = LARGE;
	minimum.psi = LARGE;
	minimum.value = LARGE;
	tempPoint.theta = 0;
	tempPoint.phi = 0;
	tempPoint.psi = 0;
	tempPoint.value = 0;
	
	current_generation_count = 1;
	reference = "";
	resetPoints();
}

void GlobalOptimisation::resetPoints()
// Emptys all vectors and makes counters zero
{
	points.empty();
	mutants.empty();
	offspring.empty();
	
	current = current_mutant = current_offspring = 0;
	bPopulationFull = bMutantsFull = bOffspringFull = false;
}

void GlobalOptimisation::setVariables(GaV ga_v)
// Sets variables for genetic algorithm run
// Inputs: string(s) - filename containing variables
{
	size = ga_v.population_size;
	generations = ga_v.generations;
	noff = ga_v.number_offspring;
	tsize = ga_v.tournament_size;
	mrate = ga_v.mutation_rate;
	
	bScreen = ga_v.print_to_screen;
	
	// Set GA variables
	//error = om.setVariables(idum, bScreen, ga_v.mutation_type);
	//if (!error) error = oo.setVariables(idum, output_content, bScreen, tsize, ga_v.mating_type, ga_v.parent_select, ga_v.fitness_type);
	om.setVariables(idum, output_content, bScreen, ga_v.mutation_type);
	oo.setVariables(idum, output_content, bScreen, tsize, ga_v.mating_type, ga_v.parent_select, ga_v.fitness_type);
	
	// MUTATION TYPES //
	if (cmpStr(ga_v.mutation_who,"offspring"))
	{
		bMChildren = true;
	}
	else if (cmpStr(ga_v.mutation_who,"parents")) 
	{
		bMChildren = false;
	}
	else
	{
#pragma omp critical
		{
			output_content->push_back("Invalid mutation option: " + ga_v.mutation_who);
		}
	}
	
	// Organise vectors ready for calculations //
	numberOfMutants = (int) round(mrate*size);
	mutants.resize(numberOfMutants);
	numberOfOffspring = (int) round(noff*size);
	offspring.resize(numberOfOffspring);	
	points.resize(size);
}

bool GlobalOptimisation::getPopulationFull()
// Check if population is full, and this if are calculations have converged
// Outputs bool if vector is full
{ 	
	if (current >= size)
	{
		// If the first time we've seen the vector full, see if converged
		if (!bPopulationFull) 
		{
			checkConverged();
		}
		bPopulationFull = true;
	}
	else 
	{
		bPopulationFull = false;
	}
	
	return bPopulationFull;
}

bool GlobalOptimisation::getMutantsFull()
// Checks if mutants is full, and then returns value
// Outputs bool if vector is full
{
	if (current_mutant >= numberOfMutants) 
	{
		bMutantsFull = true;
	}
	else 
	{
		bMutantsFull = false;
	}
	
	return bMutantsFull;
}

bool GlobalOptimisation::getOffspringFull()
// Checks if offspring is full, and then returns value
// Returns bool if vector offspring is full
{
	if (current_offspring >= numberOfOffspring) 
	{
		bOffspringFull = true;
	}
	else 
	{
		bOffspringFull = false;
	}
	
	return bOffspringFull;
}

void GlobalOptimisation::setPoint(const RotationPoint &l) 
// Sets point value in results arrays. Deals with population, mutants and offspring
// Inputs : RotationPoint(l) - value to add
{ 
	if (!bPopulationFull)
	{
		points[current] = l;
		next();
	}
	else if (!bMutantsFull) 
	{	
		mutants[current_mutant] = l;
		nextMutant();
	}
	else if (!bOffspringFull) 
	{	
		offspring[current_offspring] = l;
		nextOffspring();
	}
	
	// Output vector sizes if screen enabled //
	if (bScreen)
	{
		printSetPoints();	
	}
}

void GlobalOptimisation::setPoints(const std::vector<RotationPoint> &rp) 
// Sets point value in results arrays. Deals with population, mutants and offspring
// Inputs : vector of RotationPoint(rp) - values to append
{ 
	if (!bPopulationFull) 
	{
		points = rp; 
		current = size = rp.size();
	}
	else if (!bMutantsFull) 
	{
		mutants = rp;
		current_mutant = numberOfMutants = rp.size();
	}
	else if (!bOffspringFull) 
	{
		offspring = rp;
		current_offspring = numberOfOffspring = rp.size();
	}
	
	// Output vector sizes if screen enabled //
	if (bScreen) 
	{
		printSetPoints();
	}
}

void GlobalOptimisation::printSetPoints()
// Prints out the current number of population members
{
	// Output vector sizes if screen enabled //
	string currentS;
	string current_mutantS;
	string current_offspringS;
	NumberToString(current,currentS);
	NumberToString(current_mutant,current_mutantS);
	NumberToString(current_offspring,current_offspringS);
	
#pragma omp critical
	{
		output_content->push_back("Set point : Parents -  " + currentS + ", Mutants - " + current_mutantS + ", Offspring - " + current_offspringS);
	}
}

RotationPoint GlobalOptimisation::getNewPoint() 
// Get point from newly calculated vectors to find out function value
// Returns : RotationPoint of value
{
	if (!bPopulationFull) 
	{
		return newPoints[current];
	}
	else if (!bMutantsFull) 
	{
		return mutants[current_mutant];
	}
	else if (!bOffspringFull) 
	{
		return offspring[current_offspring];
	}
	else 
	{
#pragma omp critical
		{
			output_content->push_back("Error: All arrays are full");
		}
		// Default if there is a mistake
		RotationPoint rp = emptyRotationPoint();
		return rp;
	}
}

std::vector<RotationPoint> GlobalOptimisation::getNewPoints()
// Get points from newly calculated vectors to find out function value
// Returns : vector of RotationPoints of value
{ 
	if (!bPopulationFull) 
	{
		return newPoints; 
	}
	else if (!bMutantsFull) 
	{
		return mutants;
	}
	else if (!bOffspringFull) 
	{
		return offspring;
	}
	else
	{
#pragma omp critical
		{
			output_content->push_back("Error: All arrays are full");
		}
		// Default if there is a mistake
		vector<RotationPoint> vRP;
		return vRP;
	}
}

void GlobalOptimisation::newPopulation()
// Calculate new population as we have all the required information
{	
	if (bScreen) 
	{
#pragma omp critical
		{
			output_content->push_back("Creating New population");
		}
	}
	
	int one = 0;
	newPoints.clear();
	
	// Add populations together //		
	points.insert(points.end(),mutants.begin(),mutants.end());
	points.insert(points.end(),offspring.begin(),offspring.end());
	
	if (bScreen) 
	{
		int num_pointsI = points.size();
		string num_pointsS;
		NumberToString(num_pointsI,num_pointsS);
		
#pragma omp critical
		{
			output_content->push_back("Points size currently : " + num_pointsS);
		}
	}
	
	// Parallelized loop until population is back to required size //
#pragma omp parallel for ordered default(none) \
private(one)	
	for (int i = 0; size > i; i++)
	{
		// ELITIST //
#pragma omp critical
		{
			one = getMinimum(points,emptyRotationPoint(),output_content);
		}
		
		newPoints.push_back(points[one]);
		points.erase(points.begin()+one);	
	}
	////////////////////////	
	if (bScreen) 
	{
		int num_pointsI = newPoints.size();
		string num_pointsS;
		string thetaS;
		string phiS;
		string psiS;
		string valueS;
		NumberToString(num_pointsI,num_pointsS);
		
#pragma omp critical
		{
			output_content->push_back("");
			output_content->push_back("Points size currently : " + num_pointsS);
			output_content->push_back("Members of new population: ");
		}
		
		for (unsigned int i = 0; i < newPoints.size(); i++)
		{
			NumberToString(newPoints[i].theta,thetaS);
			NumberToString(newPoints[i].phi,phiS);
			NumberToString(newPoints[i].psi,psiS);
			NumberToString(newPoints[i].value,valueS);
			
#pragma omp critical
			{
				output_content->push_back(thetaS + tab + phiS + tab + psiS + tab + valueS);
			}
		}
		
#pragma omp critical
		{
			output_content->push_back("");
		}
	}
	
	resetPoints();
	points.resize(size); // Resize to deal with change from adding mutants and offspring to same vector
}

void GlobalOptimisation::checkConverged()
// Look to see if search has finished
// Takes no inputs or outputs but complete check
{
	int low = getMinimum(points,emptyRotationPoint(),output_content); // Current min
	
	//cout << low << endl;
	
	// First run //
	if ((minimum.theta == LARGE) &&
		(minimum.psi == LARGE) &&
		(minimum.phi == LARGE) &&
		(minimum.value == LARGE))
	{
		minimum = points[low];
	}
	else if (minimum.value > points[low].value)
		// New point better than previous one //
	{
		minimum = points[low];
		current_generation_count = 1;
	}
	else 
	{
		current_generation_count++;
	}
	
	if (bScreen)
	{
		printMinimum();
	}
	
	// See if repeats matches total generations target //
	if (current_generation_count > generations) 
	{
		bConverged = true;
	}
	else
		// Calculate mutants and offspring if not //
	{
		offspring = oo.getOffspringPoints(numberOfOffspring,getSearchLimits(),points);
		
		if (bMChildren) 
		{
			mutants = om.getMutantPoints(numberOfMutants,current_generation_count,getSearchLimits(),offspring);
		}
		else 
		{
			mutants = om.getMutantPoints(numberOfMutants,current_generation_count,getSearchLimits(),points);
		}
	}
}

LinearStruct GlobalOptimisation::minimisationPoints(const Direction &d)
// Organise linearStruct ready for running of possible rotations in local optimisation
//		 (Directions) - Search Vectors we can use
// Outputs linearStruct once organised
{	
	LinearStruct linearStruct = getSearchLimits();
	
	linearStruct.theta_max = getTempPoint().theta + (d.theta*linearStruct.theta_step);
	linearStruct.phi_max = getTempPoint().phi + (d.phi*linearStruct.phi_step);
	linearStruct.psi_max = getTempPoint().psi + (d.psi*linearStruct.psi_step);
	
	linearStruct.theta_min = getTempPoint().theta - (d.theta*linearStruct.theta_step);
	linearStruct.phi_min = getTempPoint().phi - (d.phi*linearStruct.phi_step);
	linearStruct.psi_min = getTempPoint().psi - (d.psi*linearStruct.psi_step);
	
	linearStruct.theta_range = linearStruct.theta_max - linearStruct.theta_min;
	linearStruct.phi_range = linearStruct.phi_max - linearStruct.phi_min;
	linearStruct.psi_range = linearStruct.psi_max - linearStruct.psi_min;
	
	linearStruct.theta_no_steps = (int) round(linearStruct.theta_range/linearStruct.theta_step);
	linearStruct.phi_no_steps = (int) round(linearStruct.phi_range/linearStruct.phi_step);
	linearStruct.psi_no_steps = (int) round(linearStruct.psi_range/linearStruct.psi_step);
	
	return linearStruct;
}

void GlobalOptimisation::printMinimum()
{
	string thetaS;
	string phiS;
	string psiS;
	string valueS;
	
	NumberToString(minimum.theta,thetaS);
	NumberToString(minimum.phi,phiS);
	NumberToString(minimum.psi,psiS);
	NumberToString(minimum.value,valueS);
	
#pragma omp critical
	{
		output_content->push_back("");
		output_content->push_back("Overall Minimum value is ( " + thetaS + " , " + phiS + " , " + psiS + " ) : " + valueS);
		output_content->push_back("");
	}
}

void GlobalOptimisation::printTempPoint(const string &front)
{
	string thetaS;
	string phiS;
	string psiS;
	string valueS;
	
	NumberToString(tempPoint.theta,thetaS);
	NumberToString(tempPoint.phi,phiS);
	NumberToString(tempPoint.psi,psiS);
	NumberToString(tempPoint.value,valueS);
	
#pragma omp critical
	{
		output_content->push_back(front + "( " + thetaS + " , " + phiS + " , " + psiS + " ) : " + valueS);
	}
}
