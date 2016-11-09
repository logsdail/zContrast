#include "Search.h"

/**
 Updates:
 23/08/2011
 - copyFromHistory updated to catch values outside search error and throw non
 critical error. PBCs will stop this happening.
 - same applied to addCurrentToHistory
 20/09/2011
 - edited setOptimisationStartingPoints and getComplete to pass number of steps back for Powell's method
 - This needs consolidation and tidying up. Remove asignment for tempPoint maybe in setOptimisationStartingPoints
 21/09/2011
 - Removed horrible local minimisation section from setOptimisationStartingPoints
 - Updated locallyMinimised to save results for local minimisation
 **/

using namespace std;

void Search::init(int *s, vector<string> *o) 
// Constructor
// Make sure our variables start false
{
	// Initiate random numbers in optimisation
	optimisation.init(s,o);
	// Clear all variables
	clearUp();
	// Set Booleans to false
	bOptimisation = bGlobal = bPeriodic = bLinear = bInTheLoop = bOptimising = false;
	historyCopyCounter = 0;
	powellStepCounter = 0;
	// Set outputfile
	output_content = o;
}

void Search::clearUp()
// Reset variables which are reused
{
	// ROTATIONS CLEAROUT //
	rotations.clear();
	current_results.clear();
	// CLEAR BOOLEAN OPERATORS //
	bComplete = bOptimised = false;
}

void Search::getVariables(SearchV seV, const Variables &v)
// Reads in variables for search and works from there
// Inputs(filename) - input file
//	 (v) - Variable class containing already gathered information
// Returns boolean of success or failure
{
	var = v;
	var.bRotate = true;
	
	bOptimisation = seV.bMinimise;
	bGlobal = seV.bGlobal_search;
	if (bGlobal) 
	{
		optimisation.setVariables(seV.ga_v);
	}
	bScreen = seV.bPrint_to_screen;
	
	bPeriodic = seV.bPeriodic;
	
	if (cmpStr(seV.minimise_type,"uni"))
	{
		miniDetails.setBUni();
	}
	else if (cmpStr(seV.minimise_type,"powell")) 
	{
		miniDetails.setBPowell();
	}
	else if (cmpStr(seV.minimise_type,"multi")) 
	{
		miniDetails.setBMulti();
	}
	else if (cmpStr(seV.minimise_type,"linear")) 
	{
		bLinear = true;
	}
	else if (!bGlobal)
	{
#pragma omp critical
		{
			output_content->push_back("No search method defined in this line: " + seV.minimise_type);
		}
		setComplete();
	}
	
	Create_Rotations rot;
	
	// Check limits are in correct order, then work out ranges
	optimisation.setSearchLimits(seV.steps);
	miniDetails.setSearchLimits(seV.steps);
	
	// Resize history vector to hold all values;
	vector<float> psi (seV.steps.psi_no_steps + 1, -1);
	vector< vector<float> > phi (seV.steps.phi_no_steps + 1, psi);
	vector< vector< vector<float> > > theta (seV.steps.theta_no_steps + 1, phi);
	historyNew = theta;
	//////////////
	
	if (bOptimisation || bGlobal)
	{
		optimisation.setReference(seV.comparison);
	} 
	
	if (bGlobal) 
	{
		optimisation.setRandomStartPoints();
	}
	
	if (bOptimisation) 
	{
		setOptimisationStartingPoints();
	}
	else 
	{
		// cout << "Starting Points" << endl;
		setStartingPoints(seV.steps);
	}
	
}

void Search::setOptimisationStartingPoints()
// Designates startings points for first time run of the class
{	
	
	if (!bOptimising)
	{
		// Random start point for Global Minimisation
		if (bGlobal)
		{
			optimisation.setTempPoint(optimisation.getNewPoint());
		}
		//////////////////////////// Local Minimisation
		else if (bOptimisation)	 
		{
			optimisation.setTempPoint(var.rotateX,var.rotateY,var.rotateZ,LARGE);
		}
		
		/////////////////////////////
		if (bScreen) 
		{
			optimisation.printTempPoint("Not Optimising: ");
		}
	} 
	else 
	{
		if (bScreen) 
		{
			optimisation.printTempPoint("Optimising: ");
		}
	}
	
	setStartingPoints(optimisation.minimisationPoints(miniDetails.getDirections()));
}

void Search::setStartingPoints(const LinearStruct linearStruct)
// Sets starting rotation coordinates
// Takes linearStruct file and uses all variables.
{
	Create_Rotations rot;
	//////////////////////////////////////
	if (bLinear) 
	{
		rotations = rot.createPoints(linearStruct);
		// cout << rotations.size() << endl;
	}
	else if (bOptimisation) 
	{
		rotations = rot.createPoints(miniDetails.getBMulti(),linearStruct,optimisation.getSearchLimits(),bPeriodic);
	}
	else if (bGlobal) 
	{
		rotations = optimisation.getNewPoints();
	}
	//////////////////////////////////////
	
	// Check to see if we have already calculated points - if so we can remove them from needing calculating
	copyFromHistory();
}	


Variables Search::getNext(const unsigned int current)
// Returns next variable setup for calculations
// Outputs Variable class ready to go
{	
	if ((current+1) == rotations.size())
	{
		bComplete = true;
	}
	
	RotationPoint *rp = &rotations[current];
	var.rotateX = rp->theta;
	var.rotateY = rp->phi;
	var.rotateZ = rp->psi;
	
	return var;
}

bool Search::getComplete()
// Checks to see if different structures are now completed. Important method and brutal to understand...
// Outputs boolean  - if complete
{
	int i;
	
	if ((bOptimisation) && (bComplete))
	{
		bOptimising = true;
		// GET MINIMUM STRUCTURE
		// i = optimisation.getMinimum(current_results,optimisation.getTempPoint());
		i = getMinimum(current_results,optimisation.getTempPoint(),output_content);
		
		// COMPARE TO LAST MINIMUM
		float tx, ty, tz;
		
		tx = fabs(current_results[i].theta - optimisation.getTempPoint().theta);
		ty = fabs(current_results[i].phi - optimisation.getTempPoint().phi);
		tz = fabs(current_results[i].psi - optimisation.getTempPoint().psi);
		
		// cout << powellStepCounter << endl;
		
		// CHECK IF THE CHANGES ARE ESSENTIALLY ZERO FROM THE LAST POINT
		if ((tx < ERRORBAR) && (ty < ERRORBAR) && (tz < ERRORBAR))
		{
			// CHANGE FLAGS //
			if (miniDetails.getBMulti())
			{
				locallyMinimised(i);
			}
			else if (miniDetails.getBUni())
			{
				miniDetails.nextUnivariate(powellStepCounter);
				// CHECK IF WE HAVE COMPLETED LOOP //
				if (miniDetails.getUniFinished())
				{
					locallyMinimised(i);
				}
				else 
				{
					powellStepCounter = 0;
					optimisation.setTempPoint(current_results[i]);	// SAVE MINIMUM FROM LAST SEARCH AS NOT MINIMISED YET//
					addCurrentToHistory();				// ADD HISTORY TO CATALOGUE //
					clearUp();							// CLEAR ALL CURRENT DATA //
					setOptimisationStartingPoints();	// SET NEW STARTING POINTS USING LAST OPTIMISATION //
				}
				
			}
			/////////////////////////////////
		}
		else
		{
			powellStepCounter++;
			optimisation.setTempPoint(current_results[i]);	// SAVE MINIMUM FROM LAST SEARCH AS NOT MINIMISED YET//
			addCurrentToHistory();							// ADD HISTORY TO CATALOGUE //
			clearUp();										// CLEAR ALL CURRENT DATA //
			setOptimisationStartingPoints();				// SET NEW STARTING POINTS USING LAST OPTIMISATION //
		}
	}
	else if (bGlobal && !bOptimisation)
	{
		// CHECK TO SEE IF WE HAVE ACQUIRED RESULTS YET //
		// THIS HELPS WITH INITIAL RUN THROUGH //
		if (current_results.size() != 0)
		{
			optimisation.setPoints(current_results);		// ADD ALL POINTS TO OPTIMISATION //
			addCurrentToHistory();							// APPEND NEW INFORMATION TO HISTORY //
			// USE OPTIMISED TO SIGNIFY COMPLETION //
			if (optimisation.getPopulationFull() && optimisation.getMutantsFull() && optimisation.getOffspringFull()) 
			{
				bOptimised = true;
			}
			else 
			{	
				clearUp();
				setStartingPoints(optimisation.getSearchLimits());
			}
		}
	}
	
	if (bGlobal && bOptimised)
	{
		// CLEAR ALL DATA //
		clearUp();
		// POPULTATION IS NOT FULL //
		if (bOptimisation && (!optimisation.getPopulationFull() || !optimisation.getMutantsFull() || !optimisation.getOffspringFull()))
		{
			setOptimisationStartingPoints();
		}
		// ELSE IF NOT CONVERGED RESULTS //
		else if (optimisation.getPopulationFull() && optimisation.getMutantsFull() && optimisation.getOffspringFull() && !optimisation.getConverged()) 
		{
			optimisation.newPopulation();					// CALCUATE NEW POPULATION //
			//bInTheLoop = true;								// NOTE WE ARE NOW UP AND RUNNING //
			// SET NEW START POINTS //			
			if (bOptimisation) 
			{
				setOptimisationStartingPoints();
			}
			else 
			{
				setStartingPoints(optimisation.getSearchLimits());
			}
		}
	}
	
	if (bLinear && bComplete) 
	{	
		return true; // WHILE COMPLETE AND NOT OPTIMISING (ONLY DOES ONE LOOP) //
	}
	else if (bOptimisation && bOptimised) 
	{
		return true; // OR WHILE OPTIMISING AND OPTIMISED //  
	}
	else if (bGlobal && optimisation.getConverged())  
	{
		return true; // OR FINISHED GLOBAL SEARCH //
	}
	else if (!bLinear && !bOptimisation && !bGlobal && bComplete) 
	{
		return true; // OR JUST A ONE OFF COMPARISON //
	}
	// ELSE WE HAVE NOT FINISHED //
	else 
	{
		return false;
	}
}

void Search::locallyMinimised(const int &i)
// Change flags and save result if completed minimisation
// Inputs - int i : position of minima in current_results
{
	bOptimising = false;
	bOptimised = true;
	// SAVE POINT FOR NEXT POPULATION //
	if (bGlobal) 
	{
		optimisation.setPoint(current_results[i]);	
	}
	else
		// Otherwise this is the minimum if we are not doing Global Search
	{
		optimisation.setMinimumPoint(current_results[i]);
	}
	
}

void Search::addCurrentToHistory() 
// Adds current repository to history so can be double checked if need be
{
	LinearStruct ls = optimisation.getSearchLimits();
	int current_results_size = current_results.size();
	
	for (signed int i = 0; i < current_results_size; i++)
	{
		//We need to add a check in here for values outside the search area
		//Otherwise with pbc turned off we get a segmentation fault
		if ((rotations[i].theta >= ls.theta_min) &&
			(rotations[i].theta <= ls.theta_max) &&
			(rotations[i].phi >= ls.phi_min) &&
			(rotations[i].phi <= ls.phi_max) &&
			(rotations[i].psi >= ls.psi_min) &&
			(rotations[i].psi <= ls.psi_max))
			// So this should pick up everything inside the search area or on the boundary
		{
			// These declarations need to be moved outside the for loop!
			int tx_steps = (int) round((current_results[i].theta - ls.theta_min) / ls.theta_step); 
			int ty_steps = (int) round((current_results[i].phi - ls.phi_min) / ls.phi_step);
			int tz_steps = (int) round((current_results[i].psi - ls.psi_min) / ls.psi_step);	
			
			//cout << tx_steps << " " << ty_steps << " " << tz_steps << endl;
			//cout << historyNew[tx_steps][ty_steps][tz_steps] << " " << current_results[i].value << endl;
			
			if (historyNew[tx_steps][ty_steps][tz_steps] == -1)
			{
				//cout << tx_steps << " " << ty_steps << " " << tz_steps << endl;
				
				historyNew[tx_steps][ty_steps][tz_steps] = current_results[i].value;
				
				// Addition to deal with Gimball lock
				// i.e. When x- and z- axes line up (at phi = -pi/2 or pi/2)
				if (bPeriodic)
				{
					if ((ty_steps == 0) || (ty_steps == ls.phi_no_steps)) // We are at a pole.
					{
						float theta_psi_combined = current_results[i].theta + current_results[i].psi;
						
						// Many views are equivalent here. If I rotate +ve around x, followed by the same -ve around z I will always have the same view
						// Therefore let's give this value to all the rotations to save time.
						
						for (signed int j = 0; j <= ls.theta_no_steps; j++ )
						{
							float theta_rotate = j * ls.theta_step;
							float psi_rotate = theta_psi_combined - theta_rotate;
							
							if (psi_rotate < ls.psi_min) 
							{
								while (psi_rotate < ls.psi_min)
								{
									psi_rotate += (ls.psi_range + 1);
								}
							}
							else if (psi_rotate > ls.psi_max) 
							{
								while (psi_rotate > ls.psi_max)
								{
									psi_rotate -= (ls.psi_range + 1);
								}
							}
							int r_tz_steps = (int) round(fmod(psi_rotate,ls.psi_step));
							if (r_tz_steps == 0)
							{
								tz_steps = (int) round(psi_rotate / ls.psi_step);
								
								// j is equivalent to the use of tx_steps above
								// we have recalculated tz_steps
								historyNew[j][ty_steps][tz_steps] = current_results[i].value;
							}
						}
					}
				}
				/////////////////////////////
			}
		}
		else
			// We'll print an error statement that we're working outside the search area
			// So the value won't be saved. Think about turning pbcs on.
		{
			printHistoryArrayError(rotations[i]);
		}
	}
}

void Search::copyFromHistory()
// Copy results from history so we do not need to recalculate
{
	LinearStruct ls = optimisation.getSearchLimits();
	
	for (unsigned int i = 0; i < rotations.size(); i++)
	{
		//cout << rotations[i].theta << " " << rotations[i].phi << " " << rotations[i].psi << endl;
		
		//We need to add a check in here for values outside the search area
		//Otherwise with pbc turned off we get a segmentation fault
		if ((rotations[i].theta >= ls.theta_min) &&
			(rotations[i].theta <= ls.theta_max) &&
			(rotations[i].phi >= ls.phi_min) &&
			(rotations[i].phi <= ls.phi_max) &&
			(rotations[i].psi >= ls.psi_min) &&
			(rotations[i].psi <= ls.psi_max))
			// So this should pick up everything inside the search area or on the boundary
		{
			// Move these declarations outside the loop Andrew!	
			int tx_steps = (int) round((rotations[i].theta - ls.theta_min) / ls.theta_step); 
			int ty_steps = (int) round((rotations[i].phi - ls.phi_min) / ls.phi_step);
			int	tz_steps = (int) round((rotations[i].psi - ls.psi_min) / ls.psi_step);
			
			//cout << tx_steps << " " << ty_steps << " " << tz_steps << endl;
			//cout << historyNew[tx_steps][ty_steps][tz_steps] << endl;
			
			if (historyNew[tx_steps][ty_steps][tz_steps] != -1)
			{
				//cout << tx_steps << " " << ty_steps << " " << tz_steps << endl;
				
				rotations[i].value = historyNew[tx_steps][ty_steps][tz_steps];
				current_results.push_back(rotations[i]);
				rotations.erase(rotations.begin()+i);
				i--;
				
				historyCopyCounter++;
			}
		}
		else
			// We'll print an error statement that we're working outside the search area
			// So the value won't be saved. Think about turning pbcs on.
		{
			printHistoryArrayError(rotations[i]);
		}
	}
	
	// Stops an infinite loop when all required values have been previously calculated
	if (rotations.size() == 0) 
	{
		bComplete = true;
	}
}

void Search::printHistoryArrayError(RotationPoint rp)
// Prints out error message for problems accessing History array
// Takes the current orientation as input
{
	string thetaS;
	string phiS;
	string psiS;
	
	NumberToString(rp.theta,thetaS);
	NumberToString(rp.phi,phiS);
	NumberToString(rp.psi,psiS);
	
#pragma omp critical
	{
		output_content->push_back("Error in adding current orientation to history");
		output_content->push_back("(" + thetaS + ", " + phiS + ", " + psiS + ") is outside the search area. Think about enabling PBCs as this value won't be saved.");
	}
}

