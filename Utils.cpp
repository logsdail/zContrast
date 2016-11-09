#include "Utils.h"

/**
Updates:
23/08/2011
- Removed tidy and moved it to another GlobalOptimisationUtils
**/

using namespace std;

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters) 
// Tokenize a line of text (with a stupid American spelling)
// Input: string(str) - Input string
//        vector(tokens) - Location of output
//	  string(delimiters) - Delimiters like spaces and tabs
{	
 	string::size_type lastPos = str.find_first_not_of(delimiters, 0); // Skip delimiters at beginning.
 	string::size_type pos     = str.find_first_of(delimiters, lastPos); // Find first "non-delimiter".
	
 	while (string::npos != pos || string::npos != lastPos)
 	{
 		tokens.push_back(str.substr(lastPos, pos - lastPos)); // Found a token, add it to the vector.
 		lastPos = str.find_first_not_of(delimiters, pos); // Skip delimiters.  Note the "not_of"
 		pos = str.find_first_of(delimiters, lastPos); // Find next "non-delimiter"
 	}
}

void point_on_3D_line(const vec3d& b_vec, const vec3d& p2, vec3d& new_point, double nu)
// Go to point on 3D line
// Inputs: vec3d(b_vec) - 3D start information
// 	   vec3d(p2) - End point
//	   vec3d(new_point) - Point on line with new information
//	   double(nu) - Displacement along line
{
	vec3d d_vec; // Work out displacement
	
	d_vec.x = p2.x - b_vec.x;
	d_vec.y = p2.y - b_vec.y;
	d_vec.z = p2.z - b_vec.z;
	
	new_point.x = b_vec.x + nu * d_vec.x;
	new_point.y = b_vec.y + nu * d_vec.y;
	new_point.z = b_vec.z + nu * d_vec.z;
}

void point_on_2D_line(const DataPoint2D &p1, const DataPoint2D &p2, DataPoint2D &new_point, double nu)
// Go to point on @D line
// Inputs: vec2d(p1) - 2D start information
//         vec2d(p2) - End point
//         vec2d(new_point) - Point on line with new information
//         double(nu) - Displacement along line
{
	DataPoint2D p3; // Get displacement
	p3.x = p2.x - p1.x;
	p3.y = p2.y - p1.y;
	
	new_point.x = p1.x + nu * p3.x;
	new_point.y = p1.y + nu * p3.y;
}

void find_white_space_end(const string &s, unsigned int &i)
//Does what it says on the tin
// Inputs: string(s) - String
//	   int(i) - Current int. Returns new white space locator at end
{
	if(i >= s.size()) return;
	
	bool isWhiteSpace = false;
	if(s[i] == ' ' || s[i] == '\t')
	{
		isWhiteSpace = true;
		i++;
	}
	else
	{
		isWhiteSpace = false;
		return;
	}
	
	while(isWhiteSpace)
		if(s[i] == ' ' || s[i] == '\t'){
			isWhiteSpace = true;
			i++;
		}
		else if( i >= s.size())
		{
			isWhiteSpace = false;
			return;
		}
		else
		{
			isWhiteSpace = false;
			return;
		}
}

void find_character_end(const string &s, unsigned int &i)
// Find last location of charcters
// Inputs: string(s) - String
//         int(i) - Current int. Returns new white space locator at end
{
	if(i >= s.size()) return;
	
	bool isChar = false;
	if(s[i] == ' ' || s[i] == '\t')
	{
		isChar = false;
		return;
	}
	else
	{
		isChar = true;
		i++;
	}
	
	while(isChar)
		if(s[i] == ' ' || s[i] == '\t')
		{
			isChar = false;
			return;
		}
		else if( i >= s.size())
		{
			isChar = false; //as its the end of the line
		}
		else
		{
			isChar = true;
			i++;
		}
}

bool cmpStr(const string &s1, const string &s2)
// Compares two strings and returns true or false
// Inputs: s1 - String 1
//         s2 - String 2		
{
	if (strcmp(s1.c_str(),s2.c_str()) == 0) 
	{
		return true;
	}
	else
	{
		return false;
	}
}

int randomNumber(int hi, int *idum)
// Scales random number to max possible
// Input: int (hi) - max possible
// Output: int - answer
{	
	// return range [0..hi-1]
	return int(randomNumber(idum)*hi); // implicit cast and truncation in return
}

float randomNumberF(float hi, int *idum)
// Scales random number to max possible
// Input: int (hi) - max possible
// Output: int - answer
{	
	// return range [0..hi-1]
	return float(randomNumber(idum)*hi); // implicit cast and truncation in return
}

// This has been copied from Numerical Recipes in C
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float randomNumber(int *idum)
// Get random number between 0 and 1. 
// As implemented in Numerical Recipes in C
// Inputs:  int seed value (must be negative)
// Returns: float
// int *idum;
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
	
	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
				mj=ma[ii];
				}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
					}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	// cout << "Random Number Generated: " << mj*FAC << endl;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
