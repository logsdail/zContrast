/*
 *  compare_work_link.cpp
 *  zContrast
 *
 *  Created by Andrew Logsdail on 30/01/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

#include "CompareWorkLink.h"

using namespace std;

Information createM(Variables &tempVar, int *seed, vector<string> *output_content)
// Takes variable and extracts data to run
// Inputs(tempVar) - Variables class with information
// Outputs Information class with all read information
{
	Information data = createM(tempVar.image_outfilename, tempVar.structure_filename, tempVar.g_iX, tempVar.g_iY, tempVar.g_fGridSize,
				   tempVar.g_fExponent, tempVar.g_fAlpha, tempVar.g_fScaler, tempVar.g_fNoise, tempVar.rotateX, tempVar.rotateY,
				   tempVar.rotateZ, tempVar.bFourier, tempVar.bCrossSection, tempVar.bScreen, tempVar.bRotate, tempVar.bSave,
				   tempVar.a_elementName, tempVar.a_atomicNumber, tempVar.a_atomicRadii, seed, output_content);
	
	return data;
}

Information createM(string &image_outfilename, const string &structure_filename, const float &g_iX, const float &g_iY,
					const float &g_fGridSize, const float &g_fExponent, const float &g_fAlpha, const float &g_fScaler,
					const float &g_fNoise, float const &rotateX, float const &rotateY, float const &rotateZ,
					bool const &bFourier, bool const &bCrossSection, bool const &bScreen, bool const &bRotate,
					const bool &bSave, const vector<string> &a_elementName, const vector<int> &a_atomicNumber,
					const vector<float> &a_atomicRadii, int *seed, vector<string> *output_content)
// Creates new Zcontrast images using data from input file
// Wow this is going to be a big descriptor!
// Inputs(image_outfilename) - output filename
//	 (structure_filename) - input filename
//	 (g_iX) - grid width
//	 (g_iY) - grid height. Both are integers
//	 (g_fGridSize) - space between points on grid
//	 (g_fExponent) - exponent on rutherford scattering
//	 (g_fAlpha) - Gaussian distributable variable
//	 (g_fScaler) - Scaler for Fourier transform
//	 (g_fNoise) - Maximum noise value on reading
//	 (rotateX) - rotation for X
//	 (rotateY) - rotation for Y
//	 (rotateZ) - rotation for Z
//	 (bFourier) - do we fourier?
//	 (bCrossSection) - do we cross section?
//	 (bScreen) - do we output to screen?
//	 (bRotate) - do we rotate?
//	 (bSave) - do we save files?
//	 (a_elementName) - element types
//	 (a_atomicNumber) - atomic numbers from parameters
//	 (a_atomicRadii) - atomic radii
// Outputs Information file once calculation is complete
{
	// INITIALISE VARIABLES //
	const string zcon="zcon";
	const string txt=".txt";
	const string xyz=".xyz";
	string sub = structure_filename.substr(structure_filename.size()-4);
	string temp ="";
	/////////////////////////
	
	// START STRUCTURE CALCULATION //
	if (bScreen)
	{
#pragma omp critical
		{
			output_content->push_back("");
			output_content->push_back("========================================================");
		}
	}
	
	CZcontrast zcontrast(seed,output_content);
	CContrastImage image;

	if (cmpStr(sub,xyz))
	{
		CCluster structure(structure_filename);
		
		
		// ROTATE STRUCTURE IF DESIRED //
		if (bRotate)
		{
			string rotateXS;
			string rotateYS;
			string rotateZS;
			
			NumberToString(rotateX,rotateXS);
			NumberToString(rotateY,rotateYS);
			NumberToString(rotateZ,rotateZS);
			
			if (bScreen) 
			{
#pragma omp critical
				{
					output_content->push_back("== Performing rotation: x = " + rotateXS + " , y = " + rotateYS + " , z = " +rotateZS);
				}
			}
			
			structure.rotate_x_axis(rotateX*DEG2RAD);
			structure.rotate_y_axis(rotateY*DEG2RAD);
			structure.rotate_z_axis(rotateZ*DEG2RAD);
			
			image_outfilename += "-x" + rotateXS;
			image_outfilename += "-y" + rotateYS;
			image_outfilename += "-z" + rotateZS;
		}
		
		zcontrast.setGridSize(g_fGridSize);
		zcontrast.setImageArea(g_iX,g_iY);        
		zcontrast.setCluster(structure);
		zcontrast.setConstituents(a_elementName,a_atomicNumber,a_atomicRadii);
		zcontrast.setExponent(g_fExponent);
		zcontrast.setAlpha(g_fAlpha);
		zcontrast.setNoise(g_fNoise);
		
		if (bScreen)
		{ 
#pragma omp critical
			{
				output_content->push_back("== Configuration Successful for: " + structure_filename);
				output_content->push_back("== Starting Image Simulation");
			}
		}
		
		// GENERATE IMAGE //
		zcontrast.Image();
		image = zcontrast.getImage();
		
		if (bSave)
		{
			if (bScreen) 
			{
#pragma omp critical
				{
					output_content->push_back("== Saving image: " + image_outfilename + ".zcon");
				}
			}
			// SAVE IMAGE
			image.saveImage(image_outfilename + ".zcon");
			image.saveImageMatrix(image_outfilename + ".txt");
		}
	} 
	else if (cmpStr(sub,zcon) || cmpStr(sub,txt))
	{
		if (bScreen) 
		{
#pragma omp critical
			{
				output_content->push_back("== Loading image from file : " + structure_filename);
			}
		}

		image.readImage(structure_filename);

		if (cmpStr(sub,txt))
		{
			image.saveImage(image_outfilename + ".zcon");
		}
	}
	else 
	{
#pragma omp critical
		{
			output_content->push_back("== Unrecognised file type : " + structure_filename);
		}
		return image.getInformation("File Error");
	}

	// GENERATE CROSS SECTION //
	if (bCrossSection && bSave)
	{
		int y_variable = (int) ((g_iY/2)*(1/g_fGridSize));
		if (bScreen) 
		{
#pragma omp critical
			{
				output_content->push_back("== Creating and writing cross section of image");
			}
		}
		image.saveCrossSection(image_outfilename + ".contrast", y_variable);
	}

	// CREATE FOURIER TRANSFORM //
	if (bFourier)
	{
		image.fourier(g_fScaler);
		if (bScreen) 
		{
			string g_fScalerS;
			NumberToString(g_fScaler,g_fScalerS);
			
			output_content->push_back("== Fourier Transform Successful, Scaler = " + g_fScalerS);
		}
		if (bSave) 
		{
			if (bScreen) 
			{
#pragma omp critical
				{
					output_content->push_back("== Saving image: " + image_outfilename + "-fourier.zcon");
				}
			}
			image.saveImage(image_outfilename + "-fourier.zcon");
		}
	}

	if (bScreen)
	{
#pragma omp critical
		{
			output_content->push_back("== Finished");
			output_content->push_back("========================================================");
			output_content->push_back("");
		}
	}

	Information data = image.getInformation(image_outfilename + ".xyz");

	return data;
}

Information checkGridsMatch(bool check, Search *search, Information file1, Information *file2,
							Variables tv, Variables *var, const bool fitting, int *seed, std::vector<std::string> *output_content)
// Check to make sure all the grids are of matching size for comparison
{
	if (check)
	{
		// NEED TO REGENERATE IMAGE IF THIS IS THE CASE SO SIZES MATCH
		float newGrid = printRegenerateGrid(&file1,file2,output_content);			
		
		tv.g_iX = var->g_iX = (file1.getM()*file1.getStep());
		tv.g_iY = var->g_iY = (file1.getN()*file1.getStep());
		tv.g_fGridSize = var->g_fGridSize = newGrid;
		
		if (!fitting) search->setVariables(tv);
		{
			file1 = createM(tv,seed,output_content);
		}
		// SET NEW DATA ARRAYS TO COMPARISON FILE //
	}
	
	return file1;
}

// SCALE RESULTS //
void scaleZValues(Information *file1,Information *file2, CompareV coV, std::vector<std::string> *output_content)
// Scalar method to match maximum z values for the two classes
// Inputs - pointer to file 1
//		  - pointer to file 2
//	      - pointer to where the results are stored so they can be updated
{
	if (coV.scale_intensity) 
	{	
		if (coV.print_to_screen) 
		{
#pragma omp critical
			{
				output_content->push_back("Scaling Z Values to match");
			}
		}
		
		file1->scale(file2->getMaxZ()/file1->getMaxZ());                       
		//results->setDataArrayOne(file1->getDataArray(),(int) file1->getM(),(int) file1->getN());
	}
	
}

double runLsf(Compare_Data *results, vector<string> *lsfResults, Search *search, const float x, const float y,
			  const float z, const string &filename1, const string &filename2, const int counter, const int m,
			  const int n, const int step, bool check, CompareV coV, int *seed, std::vector<std::string> *output_content)
// Works out the lsf for this matching pair, and then stores all relevant information. Lots of inherited variables
// Inputs - results - inherited information on images
//		  - lsfResults - inherited results array
//		  - search - inherited search information
//		  - x - X rotation for optimisation
//		  - y - Y rotation for optimisation
//		  - z - Z rotation for optimisation
//		  - filename1 - filename1 for saving
//		  - filename2 - filename2 for saving
//		  - counter - comparison number
//		  - m - width of grid
//		  - n - height of grid
//		  - step - step size on grid
//		  - check - if optimising and correct matching word covariance
{
	double answer = 0;
	
	if (coV.lsf)
	{
		string temp = "";
		// Get results
		answer = results->lsf(coV.save_lsf_difference);
		// Print to screen
		if (coV.print_to_screen) 
		{
			string answerS;
			NumberToString(answer,answerS);
#pragma omp critical
			{
				output_content->push_back("Least Squares: " + answerS);
			}
		}
		// Save results
		if (coV.save_results)
		{
			NumberToString(answer,temp);
			temp += tab + filename1 + tab + filename2;
#pragma omp critical
                        {
				lsfResults->push_back(temp);
			}
			// Add in method to save difference, to help visualise gaps //
			// Use sparingly as IO's will kill the program
			if (coV.save_lsf_difference)
			{
				Information difference(seed,output_content);
				string diff_filename = "output_";
				string count;
				NumberToString(counter,count);
				diff_filename += count + ".diff";
				
				difference.saveImage(diff_filename,results->getLSFDifference(),m,n,step);
			}
		}
		// Store for optimisation
		if (check) 
		{
			search->setOptimisationValue(x,y,z,answer);
		}
	}
	
	return answer;
}

double runCovariance(Compare_Data *results, vector<string> *covarianceResults, Search *search, const float x, const float y,
					 const float z, const string &filename1,const string &filename2, const float file1_mean, 
					 const float file2_mean, bool check, CompareV coV, int *seed, std::vector<std::string> *output_content)
// Works out the covariance for this matching pair, and then stores all relevant information. Lots of inherited variables
// Inputs - results - inherited information on images
//		  - covarianceResults - inherited results array
//		  - search - inherited search information
//		  - x - X rotation for optimisation
//		  - y - Y rotation for optimisation
//		  - z - Z rotation for optimisation
//		  - filename1 - filename1 for saving
//		  - filename2 - filename2 for saving
//		  - file1_mean - Mean Z value for file 1
//		  - file2_mean - Mean Z value for file 2
//		  - check - if optimising and correct matching word covariance

{
	double answer = 0;
	
	if (coV.covariance)
	{
		
		string temp = "";
		// Get result
		answer = results->covariance() - (file1_mean*file2_mean);
		// Print to screen
		if (coV.print_to_screen) 
		{
			string answerS;
			NumberToString(answer,answerS);
#pragma omp critical
			{
				output_content->push_back("Covariance: " + answerS);
			}
		}
		// Save to file
		if (coV.save_results)
		{	
			NumberToString(answer,temp);
			temp += tab + filename1 + tab + filename2;
#pragma omp critical
                        {
				covarianceResults->push_back(temp);
			}
		}
		// Store for optimisation
		if (check) 
		{
			search->setOptimisationValue(x,y,z,answer);
		}
		
	}
	
	return answer;
}

void setComparativeData(const int &counter, const string &filename1, Information *file1, const string &filename2, 
						Information *file2, Compare_Data *results, CompareV coV, std::vector<std::string> *output_content)
// Display details to screen
{
	results->setDataArrayOne(file1->getDataArray(),(int) file1->getM(),(int) file1->getN());
	results->setCount(file1->getCount());
	results->setDataArrayTwo(file2->getDataArray());
	
	if (coV.print_to_screen)
	{
		float f1step = file1->getStep();
		float f1minx = file1->getMinX();
		float f1maxx = file1->getMaxX();
		float f1miny = file1->getMinY();
		float f1maxy = file1->getMaxY();
		float f1minz = file1->getMinZ();
		float f1maxz = file1->getMaxZ();
		float f1average = file1->getMean();
		float f1m = file1->getM();
		float f1n = file1->getN();
		float f1count = file1->getCount();
		
		float f2step = file2->getStep();
		float f2minx = file2->getMinX();
		float f2maxx = file2->getMaxX();
		float f2miny = file2->getMinY();
		float f2maxy = file2->getMaxY();
		float f2minz = file2->getMinZ();
		float f2maxz = file2->getMaxZ();
		float f2average = file2->getMean();
		float f2m = file2->getM();
		float f2n = file2->getN();
		float f2count = file2->getCount();
		
		string counterS;
		string f1stepS;
		string f1minxS;
		string f1maxxS;
		string f1minyS;
		string f1maxyS;
		string f1minzS;
		string f1maxzS;
		string f1averageS;
		string f1mS;
		string f1nS;
		string f1countS;
		string f2stepS;
		string f2minxS;
		string f2maxxS;
		string f2minyS;
		string f2maxyS;
		string f2minzS;
		string f2maxzS;
		string f2averageS;
		string f2mS;
		string f2nS;
		string f2countS;
		
		NumberToString(counter,counterS);
		NumberToString(f1step,f1stepS);
		NumberToString(f1minx,f1minxS);
		NumberToString(f1maxx,f1maxxS);
		NumberToString(f1miny,f1minyS);
		NumberToString(f1maxy,f1maxyS);
		NumberToString(f1minz,f1minzS);
		NumberToString(f1maxz,f1maxzS);
		NumberToString(f1average,f1averageS);
		NumberToString(f1m,f1mS);
		NumberToString(f1n,f1nS);
		NumberToString(f1count,f1countS);
		NumberToString(f2step,f2stepS);
		NumberToString(f2minx,f2minxS);
		NumberToString(f2maxx,f2maxxS);
		NumberToString(f2miny,f2minyS);
		NumberToString(f2maxy,f2maxyS);
		NumberToString(f2minz,f2minzS);
		NumberToString(f2maxz,f2maxzS);
		NumberToString(f2average,f2averageS);
		NumberToString(f2m,f2mS);
		NumberToString(f2n,f2nS);
		NumberToString(f2count,f2countS);
		
#pragma omp critical
		{
			// PRINT DETAILS TO SCREEN (DEBUG) //
			output_content->push_back("****** COMPARISON NO: " + counterS + " ******");
			output_content->push_back("");
			output_content->push_back("File 1: " + filename1);
			output_content->push_back("Step: " + f1stepS);
			output_content->push_back("Max X: " + f1maxxS);
			output_content->push_back("Min X: " + f1minxS);
			output_content->push_back("Max Y: " + f1maxyS);
			output_content->push_back("Min Y: " + f1minyS);
			output_content->push_back("Max Z: " + f1maxzS);
			output_content->push_back("Min Z: " + f1minzS);
			output_content->push_back("Average Z: " + f1averageS);
			output_content->push_back("M: " + f1mS);
			output_content->push_back("N: " + f1nS);
			output_content->push_back("Count: " + f1countS);
			output_content->push_back("");
			output_content->push_back("File 2: " + filename2);
			output_content->push_back("Step: " + f2stepS);
			output_content->push_back("Max X: " + f2maxxS);
			output_content->push_back("Min X: " + f2minxS);
			output_content->push_back("Max Y: " + f2maxyS);
			output_content->push_back("Min Y: " + f2minyS);
			output_content->push_back("Max Z: " + f2maxzS);
			output_content->push_back("Min Z: " + f2minzS);
			output_content->push_back("Average Z: " + f2averageS);
			output_content->push_back("M: " + f2mS);
			output_content->push_back("N: " + f2nS);
			output_content->push_back("Count: " + f2countS);
			output_content->push_back("");	
			/////////////////////////////////////
		}
	}
}

float printRegenerateGrid(Information *file1, Information *file2, std::vector<std::string> *output_content)
// Error message re: grid size
{
	float newGrid = (file1->getM()*file1->getStep())/file2->getM();
	float f1m = file1->getM();
	float f1n = file1->getN();
	float f1step = file1->getStep();
	float f2m = file2->getM();
	float f2n = file2->getN();
	float f2step = file2->getStep();
	string f1mS;
	string f1nS;
	string f1stepS;
	string f2mS;
	string f2nS;
	string f2stepS;
	string newGridS;
	
	NumberToString(f1m,f1mS);
	NumberToString(f1n,f1nS);
	NumberToString(f1step,f1stepS);
	NumberToString(f2m,f2mS);
	NumberToString(f2n,f2nS);
	NumberToString(f2step,f2stepS);
	NumberToString(newGrid,newGridS);
	
#pragma omp critical
	{
		output_content->push_back("!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!");	
		output_content->push_back("Grids are not of matching size so resizing...");
		output_content->push_back("IF YOU ARE WORKING WITH DECIMAL POINTS THIS");
		output_content->push_back("MAY CAUSE PROBLEMS. RESIZE IMAGE IF POSSIBLE.");
		output_content->push_back("Recommend changing this in input file!");
		output_content->push_back("Currently:");
		output_content->push_back("File 1: Width = "+ f1mS + ", Height = " + f1nS + ", Grid Size = " + f1stepS);
		output_content->push_back("File 2: Width = " + f2mS + ", Height = " + f2nS + ", Grid Size = " + f2stepS);
		output_content->push_back("After resize:");
		output_content->push_back("File 1: Width = " + f1mS + ", Grid Size = " + newGridS);
		output_content->push_back("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		output_content->push_back("");
	}
	
	return newGrid;
}

