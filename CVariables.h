#ifndef VARIABLES_H
#define VARIABLES_H

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 02/08/2011
 - This class doesn't need output_content so it has been removed
 from the constructurs
 **/

// #include <cstdlib>
// #include <vector>
// #include "Utils.h"
#include "GlobalOptimisationUtils.h"
// #include "Structures.h"

// This is going to read in the Zcontrast variables, and have them to return.
// Saves rereading every other location
class Variables{
	
public:
	Variables() 
	{idum = NULL;}
	
	~Variables() {;}
	
	//Variables(int *s, std::vector<std::string> *o) 
	//{ init("",s,o);}
	
	//Variables(const std::string &filename, std::vector<std::string> *o) 
	//{ init(filename,idum,o); }
	
	//Variables(const std::string &filename, int *s, std::vector<std::string> *o) 
	//{ init(filename,s,o); }
	
	//void init(const std::string &filename, int *s, std::vector<std::string> *o);
	
	Variables(int *s) 
	{ init("",s);}
	
	Variables(const std::string &filename) 
	{ init(filename,idum); }
	
	Variables(const std::string &filename, int *s) 
	{ init(filename,s); }
	
	void init(const std::string &filename, int *s);
	
	void setVariables(CreateV cr_v);
	
	void setStructureFilename(const std::string &filename);
	
	std::string image_outfilename;
	std::string structure_filename;
	
	float g_iX;
	float g_iY; 
	float g_fGridSize;
	float g_fExponent;
	float g_fAlpha;
	float g_fScaler;
	float g_fNoise;
	float rotateX;
	float rotateY;
	float rotateZ; 
	
	bool bFourier; 
	bool bCrossSection; 
	bool bScreen; 
	bool bRotate; 
	bool bSave;
	bool bParameters;
	
	std::vector<std::string> a_elementName;
	
	std::vector<int> a_atomicNumber;
	
	std::vector<float> a_atomicRadii;

private:
	int *idum; // Pointer for random number generator
	// std::vector<std::string> *output_content;
	
};
#endif
