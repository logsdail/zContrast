#include "CreateRotations.h"

/**
 26/08/2011
 - Updated pbc check to remove any rotations which are outside the search area
 after some problems with testing powells method
 **/

using namespace std;

vector<RotationPoint> Create_Rotations::createPoints(const bool &bMulti, const LinearStruct &linearStruct, const LinearStruct &limits, const bool &pbc)
// Checks if we are creating univariate rotation points, or multivariate
// Inputs - bMulti - Defines if multivariate
//		  - linearStruct - Mins and Maxs for search
//		  - limits - limits for global search
// Returns vector of rotation points
{
	vector<RotationPoint> rotations;
	
	if (bMulti)
	{
		rotations = createPointsLinear(linearStruct);
	}
	else 
	{
		rotations = createPointsUnivariate(linearStruct);
	}
	
	// This will tidy up all points and make sure we don't go outside search area
	rotations = pbc_check(pbc,rotations,linearStruct,limits);
	
	return rotations;
}

vector<RotationPoint> Create_Rotations::createPoints(const LinearStruct &linearStruct)
// Create multivariate rotation points to suit needs
// Inputs - linearStruct - Mins and Maxs for search
// Returns vector of rotation points
{
	
	return createPointsLinear(linearStruct);
	
	/** 
	 The following sections have been removed due to discrepancies in point distribution
	 - triangulation
	 - circular based points
	 
	 The main problem being evenly spreading points on the surface of a sphere. A square grid wrapped around
	 a sphere does not give a uniform grid i.e.
	 
	 theta = R[0,2pi], phi = R[0, pi]
	 
	 due to proximity of points near the poles. Instead, we need to approach it with U,V = R[0,1] and then:
	 
	 theta = 2.pi.U
	 phi = cos^-1(2.v - 1)
	 
	 See http://mathworld.wolfram.com/SpherePointPicking.html
	 
	 **/
}

vector<RotationPoint> Create_Rotations::createPointsUnivariate(const LinearStruct &linearStruct)
// Create next steps along linear search for minima
// Inputs - linearStruct - Mins and Maxs for this next step. Origin is calculated as average
// Returns vector with three rotation points: Min, Max and Middle
{
	RotationPoint temp = emptyRotationPoint();
	vector<RotationPoint> rotations;
	
	// MINIMUM VALUE //
	temp.theta = linearStruct.theta_min;
	temp.phi = linearStruct.phi_min;
	temp.psi = linearStruct.psi_min;
	rotations.push_back(tidy(temp));
	
	// AVERAGE VALUE //
	temp.theta = (linearStruct.theta_min+linearStruct.theta_max)/2;
	temp.phi = (linearStruct.phi_min+linearStruct.phi_max)/2;
	temp.psi = (linearStruct.psi_min+linearStruct.psi_max)/2;
	rotations.push_back(tidy(temp));
	
	// MAXIMUM VALUE //
	temp.theta = linearStruct.theta_max;
	temp.phi = linearStruct.phi_max;
	temp.psi = linearStruct.psi_max;
	rotations.push_back(tidy(temp));
	
	return rotations;
}

vector<RotationPoint> Create_Rotations::createPointsLinear(const LinearStruct linearStruct)
// Returns next variable setup for calculations
// Inputs - linearStruct - All the variables associated with rotations
//	    x_temp - Current X location
//	    y_temp - Current Y location
//	    z_temp - Current Z location
// Outputs Variable class ready to go
{
	RotationPoint temp = emptyRotationPoint();
	vector<RotationPoint> rotations;
	
	// Define loop limits //
	int theta_max_steps = linearStruct.theta_no_steps; 
	int phi_max_steps = linearStruct.phi_no_steps;
	int psi_max_steps = linearStruct.psi_no_steps;
	
	/**
	 #pragma omp parallel for ordered default(none) \
	 private(temp,linearStruct) \
	 shared(rotations,theta_max_steps,phi_max_steps,psi_max_steps) 
	 **/
	for (signed int tx = 0; tx <= theta_max_steps; tx++ )
	{
		for (signed int ty = 0; ty <= phi_max_steps; ty++ )
		{
			for (signed int tz = 0; tz <= psi_max_steps; tz++ )
			{
				temp.theta = linearStruct.theta_min + (tx * linearStruct.theta_step);
				temp.phi = linearStruct.phi_min + (ty * linearStruct.phi_step);
				temp.psi = linearStruct.psi_min + (tz * linearStruct.psi_step);
				
				rotations.push_back(tidy(temp));
			}
		}
	}
	
	return rotations;
	
}

vector<RotationPoint> Create_Rotations::pbc_check(const bool pbc, vector<RotationPoint> rotations, const LinearStruct linearStruct, const LinearStruct &limits)
{
	int rotations_size = rotations.size();
	
    // We are going to update this to remove rotations which are outside the search area should they appear.
	// Deal with pbc issues
	if (!((linearStruct.theta_min == limits.theta_min) &&
		  (linearStruct.theta_max == limits.theta_max) &&
		  (linearStruct.phi_min == limits.phi_min) &&
		  (linearStruct.phi_max == limits.phi_max) &&
		  (linearStruct.psi_min == limits.psi_min) &&
		  (linearStruct.psi_max == limits.psi_max)))
	{
		LinearStruct l = limits;
		
		// In PBC loop.
#pragma omp parallel for ordered default(none) \
shared(rotations,rotations_size,l)
		for (signed int i = 0; i < rotations_size ; i++)
		{
			bool deleted = false;
			
			if (rotations[i].phi < l.phi_min)
			{
				if (pbc)
				{
					while (rotations[i].phi < l.phi_min)
					{
						// Check if we are using mirror planes
						if ((l.theta_range) > 179) rotations[i].theta += round((l.theta_range + 1) / 2);
						rotations[i].phi = l.phi_min + (l.phi_min - rotations[i].phi);
					}
				}
				else if (!deleted)
				{
					rotations.erase(rotations.begin()+i);
					deleted = true;
					i--;
				}
			}
			
			else if (rotations[i].phi > l.phi_max)
			{
				if (pbc)
				{
					while (rotations[i].phi > l.phi_max)
					{
						// Check if we are using mirror planes
						if ((l.theta_range) > 179) rotations[i].theta += round((l.theta_range + 1) / 2);
						rotations[i].phi = l.phi_max + (l.phi_max - rotations[i].phi);
					}
				}
				else if (!deleted)
				{
					rotations.erase(rotations.begin()+i);
					deleted = true;
					i--;
				}
			}
			
			
			if (rotations[i].theta < l.theta_min) 
			{
				if (pbc)
				{
					while (rotations[i].theta < l.theta_min)
					{
						rotations[i].theta += (l.theta_range + 1);
					}
				}
				else if (!deleted)
				{
					rotations.erase(rotations.begin()+i);
					deleted = true;
					i--;
				}
			}
			
			else if (rotations[i].theta > l.theta_max) 
			{
				if (pbc)
				{
					while (rotations[i].theta > l.theta_max)
					{
						rotations[i].theta -= (l.theta_range + 1);
					}
				}
				else if (!deleted)
				{
					rotations.erase(rotations.begin()+i);
					deleted = true;
					i--;
				}
			}
			
			
			if (rotations[i].psi < l.psi_min) 
			{
				if (pbc)
				{
					while (rotations[i].psi < l.psi_min)
					{
						rotations[i].psi += (l.psi_range + 1);
					}
				}
				else if (!deleted)
				{
					rotations.erase(rotations.begin()+i);
					deleted = true;
					i--;
				}
			}
			
			else if (rotations[i].psi > l.psi_max) 
			{
				if (pbc)
				{
					while (rotations[i].psi > l.psi_max)
					{
						rotations[i].psi -= (l.psi_range + 1);
					}
				}
				else if (!deleted)
				{
					rotations.erase(rotations.begin()+i);
					deleted = true;
					i--;
				}
			}
		}
	}
	
	return rotations;
}
