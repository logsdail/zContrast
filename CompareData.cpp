#include "CompareData.h"

/**
 25/10/2011
 - Changed the way the LSF is calculated in an attempt to reduce segmentation faults.
 **/

using namespace std;

Compare_Data::Compare_Data() 
// Constructor sets values to 0
{
	count = 0;
	m = 0;
	n = 0;
	
	// Set pointers to NULL
	output_content = NULL;
}

double Compare_Data::lsf(const bool d)
// Calculates least squares fitting total for two data sets
// To make sure this isn't a massive number we have divded by the numbers of points to give an LSF average.
// Outputs: float - result
{
	// Clear out vector for saving differences
	if (d) 
	{
		difference.clear();
	}
	
	if (count != 0)
	{
		double result = 0;
		for (int i = 0; i < count; i++)
		{
			// float diff = pow(dataOne[i]-dataTwo[i],2);
			float diff = dataOne[i]-dataTwo[i];
			diff *= diff;
			if (d)
			{
				difference.push_back(diff);
			}
			result += diff;
		}
		return (double) result/count; // True use would return just result.
	}
	else
	{
#pragma omp critical
		{
		output_content->push_back("Error: Total count is not defined");
		}
		return 0;
	}
}

double Compare_Data::covariance()
// Calculates part of covariance for two data arrays
// Covariance can be defined as (M[a*b]/N)-M[a]M[b]
// This is what we are working towards here, though subtracting M[a]M[b] needs to be done afterwards
// So actually we just give the mean of a*b
// Output: float - result
{
	if (count != 0)
	{
		double result = 0;
		for (int i = 0; i < count; i++)
		{
			result+=dataOne[i]*dataTwo[i];
		}
		return (double) result/count; // True use would return just result.
	}
	else
	{
#pragma omp critical
		{
		output_content->push_back("Error: Total count is not defined");
		}
		return 0;
	}
}

vector<float> Compare_Data::getLSFDifference()
// Returns image array of numerical differences so areas for concern can be highlighted.
{
	if (difference.size() > 0)
	{
		return difference;
	}
	else 
	{
		vector<float> d(0);
		return d;
	}

}
