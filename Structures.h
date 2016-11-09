#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>

#define ERRORBAR 0.00005
#define PI 3.141592653589793238462643383279502884197 //BIG PIE :)
#define DEG2RAD PI/180.0
#define RAD2DEG 180.0/PI
#define LARGE 99999999

/**
 24/08/2011
 - Shrunk LARGE by a factor of 100 due to errors.
 - Removed identification of lsf and cov outputs, as this is now automatically handled.
 **/

struct Direction
{
	float theta;
	float phi;
	float psi;
};

struct RotationPoint
{
	float theta;
	float phi;
	float psi;
	float value;
};

struct LinearStruct
{
	float theta_min;
	float theta_max;
	float theta_step;
	float theta_range;
	float phi_min;
	float phi_max;
	float phi_step;
	float phi_range;
	float psi_min;
	float psi_max;
	float psi_step;
	float psi_range;
	
	int theta_no_steps;
	int phi_no_steps;
	int psi_no_steps;
};

struct MainV
// Main Program Variables
{
	std::string function;
	int seed;
	int processors;
};

struct CreateV
// Creation Variables
{
	std::string xyz_filename;
	
	float width;   
	float height;          
	float gridsize;           
	float rutherford;           
	float gaussian;
	float gaussian_fit;
	float fourier_scalar;
    float noise;  
	float rotate_x;
	float rotate_y;
	float rotate_z;
	
	bool parameters_file;     
	bool fourier;            
	bool save_cross_section;   
	bool print_to_screen;   
    bool rotate;
	bool rotateRandom;
	bool save_structure;
	
	int e_jNumberOfElements;	// Number of different elements in parameter file
	
	std::vector<std::string> a_elementName;
	std::vector<int> a_atomicNumber;
	std::vector<float> a_atomicRadii;
};

struct CompareV
// Numerical Comparison Techniques
{
	std::string filename1;
	std::string filename2;
	//std::string lsf_file;
	//std::string covariance_file; 
	
	bool scale_intensity;
	bool print_to_screen;      
	bool lsf;                  
	bool covariance;              
	bool save_results;        
	bool centre_image;
	bool save_lsf_difference;
	bool translate;
	//bool translate_search;
	
	float translate_x;
	float translate_y;
};

struct GaV
// Genetic Algorithm Variables
{
	int population_size;     
	int generations;  
	int tournament_size;

	float number_offspring;     
	float mutation_rate; 
	
	std::string mutation_type;
	std::string mutation_who;    
	std::string fitness_type;        
	// std::string population_select;    
	std::string mating_type;         
	std::string parent_select; 
	
	bool print_to_screen;
};

struct SearchV
// Search Variables
{
	std::string comparison;
	std::string minimise_type;
	
	bool bMinimise;
	bool bGlobal_search;
	bool bPeriodic;
	bool bPrint_to_screen;
	
	LinearStruct steps;
	GaV ga_v;
};

#endif
