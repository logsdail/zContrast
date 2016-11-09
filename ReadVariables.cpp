/*
 *  readVariables.cpp
 *  zContrast
 *
 *  Created by Andrew Logsdail on 20/01/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

/**
 24/08/2011
 - Removed identification of lsf and cov outputs, as this is now automatically handled.
 **/

#include "ReadVariables.h"

using namespace std;

ReadVariables::ReadVariables()
// Let's initialise all the variables to their defaults here
{
	// Default variables for main program
	ma_v.function = "";
	ma_v.seed = 0;
	ma_v.processors = 1;
	// Default variables for creation
	cr_v.xyz_filename = "";
	cr_v.width = 0;   
	cr_v.height = 0;          
	cr_v.gridsize = 0;           
	cr_v.rutherford = 0;           
	cr_v.gaussian = 0;
	cr_v.gaussian_fit = 0.1;
	cr_v.fourier_scalar = 0;
    cr_v.noise = 0;                
	cr_v.parameters_file = false;
	cr_v.fourier = false;
	cr_v.save_cross_section = false;   
	cr_v.print_to_screen = false;
    cr_v.rotate = false;
	cr_v.rotateRandom = false;
	cr_v.save_structure = false;
	cr_v.rotate_x = cr_v.rotate_y = cr_v.rotate_z = 0;             
	cr_v.a_elementName.clear();
	cr_v.a_atomicNumber.clear();
	cr_v.a_atomicRadii.clear();
	// Default variables for comparing
	co_v.filename1 = "";
	co_v.filename2 = "";
	//co_v.lsf_file = "";
	//co_v.covariance_file = "";
	co_v.scale_intensity = false;
	co_v.print_to_screen = false;      
	co_v.lsf = true;                  
	co_v.covariance = false;               
	co_v.save_results = true;        
	co_v.centre_image = false;
	co_v.translate = false;
	//co_v.translate_search = false;
	co_v.save_lsf_difference = false;
	co_v.translate_x = 0;
	co_v.translate_y = 0;
	// Default variables for search
	// se_v.search_type = "";
	se_v.comparison = "";
	se_v.minimise_type = "";
	se_v.bMinimise = false;
	se_v.bGlobal_search = false;
	se_v.bPeriodic = false;
	se_v.steps.theta_max = 0;
	se_v.steps.theta_min = 0;
	se_v.steps.theta_step = 1;
	se_v.steps.phi_max = 0; // For Phi = 0, pi, only study necessary is theta = 0.
	se_v.steps.phi_min = 0;
	se_v.steps.phi_step = 1;
	se_v.steps.psi_max = 0;
	se_v.steps.psi_min = 0;
	se_v.steps.psi_step = 1;
	se_v.bPrint_to_screen = false;
	// Default variables for GA
	se_v.ga_v.population_size = 10;     
	se_v.ga_v.generations = 8;       
	se_v.ga_v.number_offspring = 0.8;     
	se_v.ga_v.tournament_size = 4;      
	se_v.ga_v.mutation_rate = 0.2;           
	se_v.ga_v.mutation_type = "static";   
	se_v.ga_v.mutation_who = "parents"; 
	se_v.ga_v.fitness_type = "exponential";        
	// se_v.ga_v.population_select = "elite";    
	se_v.ga_v.mating_type = "uniform";         
	se_v.ga_v.parent_select = "tournament";  
	se_v.ga_v.print_to_screen = false;
}

bool ReadVariables::openFile(const string &filename)
// Open pipe for reading
// Inputs(filename) - String of filename
{
	
	ifstream inData;
	inData.open(filename.c_str(),ios::in);
	bool outcome = EXIT_SUCCESS;
	
	if (!inData)
	{
		cout << "No input file for variables: " << filename << endl;
		return EXIT_FAILURE;
	}
	else
	{
		outcome = getVariables(inData);
		inData.close();
		
		return outcome;
	}
	
	return EXIT_SUCCESS;
	
}

bool ReadVariables::getVariables(ifstream &inData)
// Read in variables and away we go!
// Inputs(ifstream) - File stream with all data
{
	string temp;
	bool check = EXIT_SUCCESS;
	
	while (!inData.eof())
	{
		vector<string> words;		
		
		getline(inData,temp);
		Tokenize(temp,words,"= \t");
		
		// Percentage sign is our comment marker
		// If we have this then just read to the end of the line and move on
		if (words.size() > 0)
		{
			if (words[0][0] == '%')
			{
				// Skip line
			}
			else
				// Otherwise we have a legitimate variable
			{
				if (cmpStr(words[0].substr(0,3),"ma_"))
				{
					check = mainVariable(words[0].substr(3,words[0].length()-3),words[1]);
				}
				else if (cmpStr(words[0].substr(0,3),"cr_"))
				{
					check = createVariable(words[0].substr(3,words[0].length()-3),words[1]);
				}
				else if (cmpStr(words[0].substr(0,3),"co_"))
				{
					check = compareVariable(words[0].substr(3,words[0].length()-3),words[1]);
				}
				else if (cmpStr(words[0].substr(0,3),"se_"))
				{
					check = searchVariable(words[0].substr(3,words[0].length()-3),words[1]);
				}
				else if (cmpStr(words[0].substr(0,3),"ga_"))
				{
					check = gaVariable(words[0].substr(3,words[0].length()-3),words[1]);
				}
			}
		}
		
		if (check == EXIT_FAILURE)
		{
			cout << "Exit Failure: Error in input file." << endl;
			cout << "Responsible expressions: " << words[0] << " = " << words[1] << endl;
			return EXIT_FAILURE;
		}
		
	}
	
	// Exit check that for search the variables are in the right order numerically
	if (cmpStr(ma_v.function,"search")) 
	{
		se_v = checkOrder(se_v);
	}
	
	return EXIT_SUCCESS;
	
}

bool ReadVariables::mainVariable(const string variable, string value)
{
	if (cmpStr(variable,"function"))
	{
		ma_v.function = value;
	}
	else if (cmpStr(variable,"seed"))
	{
		StringToNumber(value,ma_v.seed);
	}
	else if (cmpStr(variable,"processors"))
	{
		StringToNumber(value,ma_v.processors);
	}
	else return EXIT_FAILURE;
	
	return EXIT_SUCCESS;
}

bool ReadVariables::createVariable(const string variable, string value)
{
	if (cmpStr(ma_v.function,"create") ||
		cmpStr(ma_v.function,"both") ||
		cmpStr(ma_v.function,"search") ||
		cmpStr(ma_v.function.substr(0,4),"fit_"))
	{
		// STRINGS //
		if (cmpStr(variable,"xyz_filename"))
		{
			cr_v.xyz_filename = value;
		}
		else if (cmpStr(variable,"parameters_file"))
		{
			cr_v.parameters_file = cmpStr(value,"yes");
			if (cr_v.parameters_file) 
			{
				return readParameters();
			}
		}
		// BOOLEANS //
		else if (cmpStr(variable,"fourier"))
		{
			cr_v.fourier = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"save_cross_section"))
		{
			cr_v.save_cross_section = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"print_to_screen"))
		{
			cr_v.print_to_screen = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"rotate_random"))
		{
			cr_v.rotateRandom = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"save_structure"))
		{
			cr_v.save_structure = cmpStr(value,"yes");
		}
		// NUMBERS //
		else if (cmpStr(variable,"rotate_x"))
		{
			StringToNumber(value,cr_v.rotate_x);
			if (cr_v.rotate_x != 0) 
			{
				cr_v.rotate = true;
			}
		}
		else if (cmpStr(variable,"rotate_y"))
		{
			StringToNumber(value,cr_v.rotate_y);
			if (cr_v.rotate_y != 0) 
			{
				cr_v.rotate = true;
			}
		}
		else if (cmpStr(variable,"rotate_z"))
		{
			StringToNumber(value,cr_v.rotate_z);
			if (cr_v.rotate_z != 0) 
			{
				cr_v.rotate = true;
			}
		}
		else if (cmpStr(variable,"width"))
		{
			StringToNumber(value,cr_v.width);
		}
		else if (cmpStr(variable,"height"))
		{
			StringToNumber(value,cr_v.height);
		}
		else if (cmpStr(variable,"gridsize"))
		{
			StringToNumber(value,cr_v.gridsize);
		}
		else if (cmpStr(variable,"rutherford"))
		{
			StringToNumber(value,cr_v.rutherford);
		}
		else if (cmpStr(variable,"gaussian"))
		{
			StringToNumber(value,cr_v.gaussian);
		}
		else if (cmpStr(variable,"gaussian_fit"))
		{
			StringToNumber(value,cr_v.gaussian_fit);
		}
		else if (cmpStr(variable,"noise"))
		{
			StringToNumber(value,cr_v.noise);
		}
		else if (cmpStr(variable,"fourier_scalar"))
		{
			StringToNumber(value,cr_v.fourier_scalar);
		}
		else 
		{
			return EXIT_FAILURE;
		}
	}
	else 
	{
		cout << "Ignoring expressions: cr_" << variable << " = " << value << " as not being used by program" << endl;
	}
	
	
	return EXIT_SUCCESS;
}

bool ReadVariables::compareVariable(const string variable, string value)
{
	if (cmpStr(ma_v.function,"compare") ||
		cmpStr(ma_v.function,"both") ||
		cmpStr(ma_v.function,"search") ||
		cmpStr(ma_v.function.substr(0,4),"fit_"))
	{
		// STRINGS //
		if (cmpStr(variable,"filename1"))
		{
			co_v.filename1 = value;
		}
		else if (cmpStr(variable,"filename2"))
		{
			co_v.filename2 = value;
		}
		// BOOLEANS //
		else if (cmpStr(variable,"scale_intensity"))
		{
			co_v.scale_intensity = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"print_to_screen"))
		{
			co_v.print_to_screen = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"lsf"))
		{
			co_v.lsf = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"covariance"))
		{
			co_v.covariance = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"save_results"))
		{
			co_v.save_results = cmpStr(value,"yes");
		}
		//else if (cmpStr(variable,"translate_search"))
		//	co_v.translate_search = cmpStr(value,"yes");
		else if (cmpStr(variable,"save_lsf_difference"))
		{
			co_v.save_lsf_difference = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"centre_image"))
		{
			co_v.centre_image = cmpStr(value,"yes");
			if (co_v.centre_image) 
			{
				co_v.translate = true;
			}
		}
		// NUMBERS //
		else if (cmpStr(variable,"translate_x"))
		{
			StringToNumber(value,co_v.translate_x);
			if (co_v.translate_x != 0)
			{
				co_v.translate = true;
			}
		}
		else if (cmpStr(variable,"translate_y"))
		{
			StringToNumber(value,co_v.translate_y);
			if (co_v.translate_y != 0) 
			{
				co_v.translate = true;
			}
		}
		else 
		{
			return EXIT_FAILURE;
		}
	}
	else
	{
		cout << "Ignoring expressions: co_" << variable << " = " << value << " as not being used by program" << endl;
	}
	
	return EXIT_SUCCESS;
}

bool ReadVariables::searchVariable(const string variable, string value)
{	
	if (cmpStr(ma_v.function,"search"))
	{
		// STRINGS //
		if (cmpStr(variable,"minimise_type"))
		{
			se_v.minimise_type = value;
		}
		else if (cmpStr(variable,"comparison"))
		{
			se_v.comparison = value;
		}
		// BOOLEANS //
		else if (cmpStr(variable,"minimise"))
		{
			se_v.bMinimise = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"global_search"))
		{
			se_v.bGlobal_search = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"print_to_screen"))
		{
			se_v.bPrint_to_screen = cmpStr(value,"yes");
		}
		else if (cmpStr(variable,"periodic"))
		{
			se_v.bPeriodic = cmpStr(value,"yes");
		}
		// NUMBERS //
		else if (cmpStr(variable,"theta_max"))
		{
			StringToNumber(value,se_v.steps.theta_max);
		}
		else if (cmpStr(variable,"theta_min"))
		{
			StringToNumber(value,se_v.steps.theta_min);
		}
		else if (cmpStr(variable,"theta_step"))
		{
			StringToNumber(value,se_v.steps.theta_step);
		}
		else if (cmpStr(variable,"phi_max"))
		{
			StringToNumber(value,se_v.steps.phi_max);
		}
		else if (cmpStr(variable,"phi_min"))
		{
			StringToNumber(value,se_v.steps.phi_min);
		}
		else if (cmpStr(variable,"phi_step"))
		{
			StringToNumber(value,se_v.steps.phi_step);
		}
		else if (cmpStr(variable,"psi_max"))
		{
			StringToNumber(value,se_v.steps.psi_max);
		}
		else if (cmpStr(variable,"psi_min"))
		{
			StringToNumber(value,se_v.steps.psi_min);
		}
		else if (cmpStr(variable,"psi_step"))
		{
			StringToNumber(value,se_v.steps.psi_step);
		}
		else 
		{
			return EXIT_FAILURE;
		}
		
	}
	else
	{
		cout << "Ignoring expressions: se_" << variable << " = " << value << " as not being used by program" << endl;
	}
	
	return EXIT_SUCCESS;
}

SearchV ReadVariables::checkOrder(SearchV seV)
// Quick method to check limits are in numerical order
// Takes Search Variables as input, returns them reordered
{
	float swap = 0;
	// Quick check to make sure limits are in the right order
	// x
	if (seV.steps.theta_min > seV.steps.theta_max)
	{
		cout << "WARNING: X MIN is greater than X MAX" << endl;
		swap = seV.steps.theta_max;
		seV.steps.theta_max = seV.steps.theta_min;
		seV.steps.theta_min = swap;
	}
	// y
	if (seV.steps.phi_min > seV.steps.phi_max)
	{
		cout << "WARNING: Y MIN is greater than Y MAX" << endl;
		swap = seV.steps.phi_max;
		seV.steps.phi_max = seV.steps.phi_min;
		seV.steps.phi_min = swap;
	}
	// z
	if (seV.steps.psi_min > seV.steps.psi_max)
	{
		cout << "WARNING: Z MIN is greater than Z MAX" << endl;
		swap = seV.steps.psi_max;
		seV.steps.psi_max = seV.steps.psi_min;
		seV.steps.psi_min = swap;
	}
	
	// Now let's work out ranges
	
	seV.steps.theta_range = seV.steps.theta_max - seV.steps.theta_min;
	seV.steps.phi_range = seV.steps.phi_max - seV.steps.phi_min;
	seV.steps.psi_range = seV.steps.psi_max - seV.steps.psi_min;
	
	seV.steps.theta_no_steps = (int) round(seV.steps.theta_range/seV.steps.theta_step);
	seV.steps.phi_no_steps = (int) round(seV.steps.phi_range/seV.steps.phi_step);
	seV.steps.psi_no_steps = (int) round(seV.steps.psi_range/seV.steps.psi_step);
	
	return seV;
}


bool ReadVariables::gaVariable(const string variable, string value)
{
	if (cmpStr(ma_v.function,"search"))
	{
		// STRINGS //
		if (cmpStr(variable,"print_to_screen"))
		{
			se_v.ga_v.print_to_screen = cmpStr(value,"yes");
		}
		// NUMBERS //
		else if (cmpStr(variable,"population_size"))
		{
			StringToNumber(value,se_v.ga_v.population_size);
		}
		else if (cmpStr(variable,"generations"))
		{
			StringToNumber(value,se_v.ga_v.generations);
		}
		else if (cmpStr(variable,"number_offspring"))
		{
			StringToNumber(value,se_v.ga_v.number_offspring);
		}
		else if (cmpStr(variable,"tournament_size"))
		{
			StringToNumber(value,se_v.ga_v.tournament_size);
		}
		else if (cmpStr(variable,"mutation_rate"))
		{
			StringToNumber(value,se_v.ga_v.mutation_rate);
		}
		else if (cmpStr(variable,"mutation_type"))
		{
			se_v.ga_v.mutation_type = value;
		}
		else if (cmpStr(variable,"mutation_who"))
		{
			se_v.ga_v.mutation_who = value;
		}
		else if (cmpStr(variable,"fitness_type"))
		{
			se_v.ga_v.fitness_type = value;
		}
		// else if (cmpStr(variable,"population_select"))
		//{
		// 	se_v.ga_v.population_select = value;
		//}
		else if (cmpStr(variable,"mating_type"))
		{
			se_v.ga_v.mating_type = value;
		}
		else if (cmpStr(variable,"parent_select"))
		{
			se_v.ga_v.parent_select = value;
		}
		else 
		{
			return EXIT_FAILURE;
		}
	}
	else
	{
		cout << "Ignoring expressions: ga_" << variable << " = " << value << " as not being used by program" << endl;
	}
	
	return EXIT_SUCCESS;
}

bool ReadVariables::readParameters()
{
	string elements_filename = "parameters.in";
	ifstream elements_file;
	elements_file.open(elements_filename.c_str(),ios::in);
	
	if(!elements_file.is_open())
	{
		cout << "Unable to open file:" << elements_filename << std::endl;
		cout << "Parameters.in name must contain in the correct order:" << endl;
		cout << "(int) Number of Elements defined in parameter file" << endl;
		cout << "(str) Element" << endl;
		cout << "(int) Atomic Number" << endl;
		cout << "(float) Atomic Radii" << endl;
		cout << "(Last three repeated for however many elements documented)" << endl;
		return EXIT_FAILURE;
	}
	
	elements_file >> cr_v.e_jNumberOfElements;
	
	cr_v.a_elementName.resize(cr_v.e_jNumberOfElements);
	cr_v.a_atomicNumber.resize(cr_v.e_jNumberOfElements);
	cr_v.a_atomicRadii.resize(cr_v.e_jNumberOfElements);
	
	int counter = 0;
	while (counter != cr_v.e_jNumberOfElements)
	{
		elements_file >> cr_v.a_elementName[counter];
		elements_file >> cr_v.a_atomicNumber[counter];
		elements_file >> cr_v.a_atomicRadii[counter];
		counter++;
	}
	
	elements_file.close();
	return EXIT_SUCCESS;
}		


