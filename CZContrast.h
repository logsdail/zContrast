#ifndef __CZCONTRAST__
#define __CZCONTRAST__

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 **/

// #include <vector>
// #include <string>
#include "CCluster.h"
#include "Utils.h"
// #include <cstdlib>
// #include "math.h"
#include "CContrastImage.h"

struct data
{
	double x;
	double y;
};

/**
 * CZcontrast()
 * Purpose: To accept a cluster object and return Image of that cluster using Z-Contrast
 *          For a given area and resolution. defined in points per angstrom
 * 
 */

class CZcontrast
{
public:
	CZcontrast() 
	{	
		// Set pointers to null
		m_pImage = NULL;
		m_pCluster = NULL;
		idum = NULL;
	}
	
	CZcontrast(int *s, std::vector<std::string> *o) 
	{ init(s,o); }
	
	~CZcontrast();
	
	void init(int *s, std::vector<std::string> *o);
	
	CContrastImage getImage()
	{return (*m_pImage); };//returns a CContrastImage object
	
	void Image();
	
	void setCluster(CCluster &cluster)
	{m_pCluster = &cluster; }; //was careless I made a copy here and then assinged the address of copy so seg fault
	
	void setImageArea(float x, float y);                     //float x and y in Angstroms
	
	void setGridSize(float grid)
	{m_fGridSize = grid;   };//1.0 = 1 angstron per grid	
	
	void setProbeSize(float size)
	{m_fProbeSize = size; };
	
	void setAlpha(float a)
	{m_fAlpha = a;};
	
	void setExponent(float a)
	{m_fExponent = a;};
	
	void setNoise(float a)
	{m_fNoise = a;};
	
	void setConstituents(std::vector<std::string> a, std::vector<int> b, std::vector<float> c)
	{m_fElementName = a; m_fAtomicNumber = b; m_fAtomicRadii = c;};

	// THIS METHOD IS NOT CURRENTLY USED //
	// void CrossSection(float x1, float y1, float x2, float y2);
	// void writeCrossSection(std::string filename);
	///////////////////////////////////////
	
private:
	CContrastImage *m_pImage;
	
	CCluster *m_pCluster; //ptr to cluster object
	
	int m_iNumX; //number of data points along the x axis
	int m_iNumY; //number of data points along the y axis
	int *idum; // Pointer for random number generator
	std::vector<std::string> *output_content;
	
	float m_fPointsPerLength; //floating value that represent points per unit length defined
	float m_fGridSize; // 1.0 = 1 angstrom per grid
	float m_fProbeSize;
	float m_fAlpha; // For the gaussian probe
	float m_fExponent; // Exponent for gaussian
	float m_fNoise; // Maximum noise levels
	
	std::vector<std::string> m_fElementName; // For consituents (i.e. parameters)
	std::vector<int> m_fAtomicNumber;
	std::vector<float> m_fAtomicRadii;
	std::vector<data> m_vData;
	
	void SearchForAtoms(std::vector<atom> vec, atom &probe, std::vector<double> &atom_count);
	
	double Gaussian(const double &distance, const int &i);
	
	double DetectIntensity(std::vector<double> atom_count);
	
	double addNoise(const double &a);
	
};

#endif
