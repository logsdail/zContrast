#include "CVariables.h"

/**
 03/06/2011
		- Changed random point to be calculated using GA Utils.
 **/

using namespace std;

void Variables::init(const string &filename, int *s)
// void Variables::init(const string &filename, int *s, std::vector<std::string> *o)
// Constructor, just calls getVariables()
// Inputs(filename) - String of filename
//		          s - Seed for random number
{
	// Initiate Random Numbers
	idum = s;
	// Set output file
	// output_content = o;
	//////////////////////////
	if (filename.size() > 0) 
	{
		setStructureFilename(filename);
	}
}

void Variables::setVariables(CreateV cr_v)
/** Set Variables for comparison.
 inputs - create variables
 **/
{
	setStructureFilename(cr_v.xyz_filename);
	g_iX = cr_v.width;
	g_iY = cr_v.height; 
	g_fGridSize = cr_v.gridsize;
	g_fExponent = cr_v.rutherford;
	g_fAlpha = cr_v.gaussian;
	g_fScaler = cr_v.fourier_scalar;
	g_fNoise = cr_v.noise;
	bFourier = cr_v.fourier; 
	bCrossSection = cr_v.save_cross_section; 
	bScreen = cr_v.print_to_screen; 
	bRotate = cr_v.rotate; 
	bSave = cr_v.save_structure;
	bParameters = cr_v.parameters_file;
	a_elementName = cr_v.a_elementName;
	a_atomicNumber = cr_v.a_atomicNumber;
	a_atomicRadii = cr_v.a_atomicRadii;
	
	// Define rotations for random variables
// #define half_rot 90
	
	if (cr_v.rotateRandom)
	{
		
		RotationPoint rp;
		LinearStruct ls;
		
		ls.theta_min = 0;
		ls.theta_max = 360;
		ls.theta_step = 1;
		ls.theta_range = ls.theta_max - ls.theta_min;
		ls.phi_min = -90;
		ls.phi_max = 90;
		ls.phi_step = 1;
		ls.phi_range = ls.phi_max - ls.phi_min;
		ls.psi_min = 0;
		ls.psi_max = 360;
		ls.psi_step = 1;
		ls.psi_range = ls.psi_max - ls.psi_min;
		
		ls.theta_no_steps = (int) round(ls.theta_range/ls.theta_step);
		ls.phi_no_steps = (int) round(ls.phi_range/ls.phi_step);
		ls.psi_no_steps = (int) round(ls.psi_range/ls.psi_step);
		
		rp = getRandomPoint(idum,ls);
		
		rotateX = rp.theta;
		rotateY = rp.phi;
		rotateZ = rp.psi;
		
		bRotate = true;
	}
	else
	{
		rotateX = cr_v.rotate_x;
		rotateY = cr_v.rotate_y;
		rotateZ = cr_v.rotate_z;
	}
}

void Variables::setStructureFilename(const string &filename)
// Sets both structure and output filenames
// Input(filename) - Structure input filename
{
	structure_filename = filename;
	image_outfilename = structure_filename.substr(0,structure_filename.size()-4);
}


