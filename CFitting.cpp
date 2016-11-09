/*
 *  fitting.cpp
 *  zContrast
 *
 *  Created by Andrew Logsdail on 26/01/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 
 02/10/2011
 - This is not parallelised. Consider in the future?
 */

#include "CFitting.h"

using namespace std;

bool Fitting::gaussian(int *idum, const CreateV &crV, const CompareV &co, vector<string> *output_content)
/** Gaussian Fitting to other image e.g. Multislice
 inputs - int pointer for random number
		- variables for creating
		- variables for comparing
 **/
{	
	// SET UP AS WOULD BE DONE IN THE CONSTRUCTOR //
	CompareV coV = co; // Compare Variables
	
	// Edits to parameters so we don't get conflicts
	coV.save_results = false;
	coV.save_lsf_difference = false;
	// SET UP COMPLETE //
	
	// Set up all
	int counter = 0;
	float gaussian_change = crV.gaussian_fit;
	double answer = 0;
	double lower = 0;
	double middle = 0;
	double  higher = 0;
	
	// LOAD IN CLASSES //
	// Work w;
	Information one(idum,output_content);
	Information two(idum,output_content);
	Compare_Data results(output_content);
	Search se;
	
	// SETUP FOR CREATE //
	Variables var(idum);
	var.setVariables(crV);
	
	// LOAD IN DATA //
	two.readFile(coV.filename2,coV);
	
	//////////////////////
	string alphaS;
	NumberToString(var.g_fAlpha,alphaS);
	
	output_content->push_back("Alpha starting value : " + alphaS);
	//////////////////////
	
	bool lower_value = true; // Check if converged
	// Run loop until we are at a minimum
	while (lower_value)
	{
		// Check how many loops we've done
		counter++;
		printConvergenceError(counter,output_content);
		
		// Set values to 0
		lower = higher = middle = 0;
		
/**		// Original working OMP constructs on first line
		#pragma omp parallel for default(none) \
		private(one, answer, se) \
		shared(var, two, idum, coV, counter, results, lower_value) \
		shared(middle, higher, lower, gaussian_change)
**/
		// Original working OMP constructs
		for (int i = 0; i < 3; i++)
		{
			answer = 0;
			Variables tv = var;

			//tv.setStructureFilename(coV.filename1);
			
			if (i == 1) 
			{
				tv.g_fAlpha += gaussian_change;
			}
			else if (i == 2) 
			{
				tv.g_fAlpha -= gaussian_change;
			}
			
			//CREATE STRUCTURE //
			one = createM(tv,idum,output_content);
			
			// CHECK GRIDS MATCH IN SIZE //
			one = checkGridsMatch((one.getCount() != two.getCount()), &se, one, &two, tv, &var, true, idum, output_content);
			
			// SCALE RESULTS //
			scaleZValues(&one,&two,coV,output_content);
			
			// SET RESULTS VARIABLES //
			setComparativeData(counter,coV.filename1,&one,coV.filename2,&two,&results,coV,output_content);
			/////////////////////////////////////

			// Temp values used to meet requirements of inputs for LSF/Covariance //
			vector<string> t;
			////////
			
			if (coV.lsf)
				// LSF //
			{
				answer = runLsf(&results,&t,&se,0,0,0,coV.filename1,coV.filename2,0,(int) one.getM(),
								(int) one.getN(), (int) one.getStep(),false,coV,idum,output_content);
			}
			else if (coV.covariance)
				// COVARIANCE //
			{
				answer = runCovariance(&results,&t,&se,0,0,0,coV.filename1,coV.filename2,one.getMean(),two.getMean(),false,coV,idum,output_content);
			}
			else
			{
				output_content->push_back("We need to use covariance or lsf for fitting");
				output_content->push_back("Please specify one or the other. Exiting fit");
				lower_value = false;
			}
			
			if (i == 0) 
			{
				middle = answer;
			}
			else if (i == 1)
			{
				higher = answer;
			}
			else if (i == 2) 
			{
				lower = answer;
			}
		}
		
		// Check if results are better than previous
		if ((middle > higher) || (middle > lower))
		{
			// Check if lower is better
			if (lower > higher)
			{
				var.g_fAlpha += gaussian_change;
				
				if (coV.print_to_screen)
				{
					NumberToString(var.g_fAlpha,alphaS);
					
					output_content->push_back("");
					output_content->push_back("Increasing Alpha by one step : " + alphaS);
					output_content->push_back("");
				}
			}
			else
			{
				var.g_fAlpha -= gaussian_change;
				
				if (coV.print_to_screen)
				{
					NumberToString(var.g_fAlpha,alphaS);
					
					output_content->push_back("");
					output_content->push_back("Decreasing Alpha by one step : " + alphaS);
					output_content->push_back("");
				}
			}
		}
		else 
			// We have reached the minimum
		{
			NumberToString(var.g_fAlpha,alphaS);
			lower_value = false;
			
			output_content->push_back("Alpha is at a minimum : " + alphaS);
		}
		
	}
	
	return EXIT_SUCCESS;
}

bool Fitting::scale(int *idum, const CreateV &crV, const CompareV &co, std::vector<std::string> *output_content)
/** Scale Fitting to other image e.g. Experimental
 inputs - int pointer for random number
		- variables for creating
		- variables for comparing
 **/
{	
	// SET UP AS WOULD BE DONE IN THE CONSTRUCTOR //
	CompareV coV = co; // Compare Variables
	
	// Edits to parameters so we don't get conflicts
	coV.save_results = false;
	coV.save_lsf_difference = false;
	// SET UP COMPLETE //
	
	// Set up all
	int counter = 0;
	double answer = 0;
	double lower = 0;
	double middle = 0;
	double higher = 0;
	
	// LOAD IN CLASSES //
	// Work w;
	Information one(idum,output_content);
	Information two(idum,output_content);
	Compare_Data results(output_content);
	Search se;
	
	// SETUP FOR CREATE //
	Variables var(idum);
	var.setVariables(crV);
	
	// LOAD IN DATA //
	two.readFile(coV.filename2,coV);
	
	//////////////////////
	int one_percent_x = (int) ceil(var.g_iX / 100);
	int one_percent_y = (int) ceil(var.g_iY / 100);
	
	string g_iXS;
	string g_iYS;
	NumberToString(var.g_iX,g_iXS);
	NumberToString(var.g_iY,g_iYS);
	
	output_content->push_back("Grid X and Y Dimension starting value : " + g_iXS + " " + g_iYS);
	//////////////////////
	
	bool lower_value = true; // Check if converged
	// Run loop until we are at a minimum
	while (lower_value)
	{
		// Check how many loops we've done
		counter++;
		printConvergenceError(counter,output_content);
		
		// Set values to 0
		lower = 0;
		middle = 0;
		higher = 0;
/**		
		#pragma omp parallel for default(none) \
		private(answer, one, results) \
		shared(one_percent_x, one_percent_y, middle, higher, lower, idum) \
		shared(var, two, coV, se, counter, lower_value)
**/		
		for (int i = 0; i < 3; i++)
		{
			answer = 0;

			Variables tv = var;
			//tv.setStructureFilename(coV.filename1);
			
			if (i == 1) 
			{
				tv.g_iX += one_percent_x;
				tv.g_iY += one_percent_y;
			}
			else if (i == 2) 
			{
				tv.g_iX -= one_percent_x;
				tv.g_iY -= one_percent_y;
			}			
			
			//CREATE STRUCTURE //
			one = createM(tv,idum,output_content);
			
			// CHECK GRIDS MATCH IN SIZE //
			one = checkGridsMatch((one.getCount() != two.getCount()), &se, one, &two, tv, &var, true, idum, output_content);
			
			// SCALE RESULTS //
			scaleZValues(&one,&two,coV,output_content);
			
			// SET RESULTS VARIABLES //
			setComparativeData(counter,coV.filename1,&one,coV.filename2,&two,&results,coV,output_content);
			/////////////////////////////////////

			// Temp values used to meet requirements of inputs for LSF/Covariance //
			vector<string> t;
			////////
			
			if (coV.lsf)
			{
				// LSF // 
				answer = runLsf(&results,&t,&se,0,0,0,coV.filename1,coV.filename2,0,(int) one.getM(),
								(int) one.getN(), (int) one.getStep(),false,coV,idum,output_content);
			}
			
			else if (coV.covariance)
			{
				// COVARIANCE //
				answer = runCovariance(&results,&t,&se,0,0,0,coV.filename1,
									   coV.filename2,one.getMean(),two.getMean(),false,coV,idum,output_content);
			}
			else
			{
				output_content->push_back("We need to use covariance or lsf for fitting");
				output_content->push_back("Please specify one or the other. Exiting fit");
				lower_value = false;
			}
			
			if (i == 0) 
			{
				middle = answer;
			}
			else if (i == 1)
			{
				higher = answer;
			}
			else if (i == 2) 
			{
				lower = answer;
			}
		}
		
		// Check if results are better than previous
		if ((middle > higher) || (middle > lower))
		{
			// Check if lower is better
			if (lower > higher)
			{
				var.g_iX += one_percent_x;
				var.g_iY += one_percent_y;
				
				if (coV.print_to_screen)
				{
					NumberToString(var.g_iX,g_iXS);
					NumberToString(var.g_iY,g_iYS);
					
					output_content->push_back("");
					output_content->push_back("Increasing X and Y dimensions by one percent : " + g_iXS + " " + g_iYS);
					output_content->push_back("");
				}
			}
			else
			{
				var.g_iX -= one_percent_x;
				var.g_iY -= one_percent_y;
				
				if (coV.print_to_screen)
				{
					NumberToString(var.g_iX,g_iXS);
					NumberToString(var.g_iY,g_iYS);
					
					output_content->push_back("");
					output_content->push_back("Decreasing X and Y dimensions by one percent : " + g_iXS + " " + g_iYS);
					output_content->push_back("");
				}
			}
		}
		else 
			// We have reached the minimum
		{
			NumberToString(var.g_iX,g_iXS);
			NumberToString(var.g_iY,g_iYS);
			
			lower_value = false;
			output_content->push_back("X and Y Dimensions are at a minimum : " + g_iXS + " " + g_iYS);
		}
		
	}
	
	return EXIT_SUCCESS;
}

void Fitting::printConvergenceError(const int &c, vector<string> *output_content)
// Gives us an error message if we are looping many times.
// Precursor for putting a better initial guess in
{
	int limit = 30;
	
	if (c > limit)
	{
		string cS;
		NumberToString(c,cS);
		
		output_content->push_back("");
		output_content->push_back("Looks like you are having convergence issues - loop " + cS);
		output_content->push_back("Make sure you haven't rotated the structure by mistake");
		output_content->push_back("");
	}	
}

/**
bool Fitting::translate(int *idum, const CreateV &crV, const CompareV &co)
{	
	// SET UP AS WOULD BE DONE IN THE CONSTRUCTOR //
	//idum = s; // Random Number seed
	//crV = cr; // Create Variables
	CompareV coV = co; // Compare Variables
	
	// Edits to parameters so we don't get conflicts
	coV.save_results = false;
	coV.save_lsf_difference = false;
	// SET UP COMPLETE //
	
	// Set up all
	int counter = 0;
	int one_percent_x, one_percent_y, total_x, total_y;
	double answer, left, right, middle;
	
	// LOAD IN CLASSES //
	// Work w;
	Information one(idum), two(idum), original(idum);
	Compare_Data results;
	Search se;
	
	// SETUP FOR CREATE //
	Variables var(idum);
	var.setVariables(crV);
	
	// LOAD IN DATA //
	original.readFile(coV.filename2,coV);
	
	//////////////////////
	one_percent_x = (int) ceil(original.getM() / 100);
	one_percent_y = (int) ceil(original.getN() / 100);
	//cout << original.getM() << " " << original.getN() << endl;
	total_x = total_y = 0;
	//////////////////////
	
	// Create input file
	one = createM(var,idum);
	//cout << one.getCount() << endl;
	
	// CHECK GRIDS MATCH IN SIZE //
	one = checkGridsMatch((one.getCount() != original.getCount()), &se, &results, one, &original, var, &var, true, idum);
	
	// SCALE RESULTS //
	scaleZValues(&one,&original,coV);
	
	bool vertical = true; // Check if converged
	bool lateral = true; //
	// Run loop until we are at a minimum
	while (vertical || lateral)
	{
		// Check how many loops we've done
		counter++;
		printConvergenceError(counter);
		
		// Set values to 0
		left = right = middle = 0;
		
		if (lateral)
		{
			
			// LATERAL SEARCH //
			for (int i = 0; i < 3; i++)
			{
				coV.translate_x = co.translate_x + total_x + ((i-1)*one_percent_x);
				two.readFile(coV.filename2,coV);
				//two.translate(total_x + ((i-1)*one_percent_x), total_y);
				
				// SET RESULTS VARIABLES //
				setComparativeData(counter,coV.filename1,&one,coV.filename2,&two,&results,coV);
				/////////////////////////////////////
				
				//cout << one.getCount() << endl;
				
				if (coV.lsf)
					// LSF // 
					answer = runLsf(&results,NULL,NULL,NULL,NULL,NULL,coV.filename1,coV.filename2,NULL,(int) one.getM(),
									(int) one.getN(),one.getStep(),false,coV,idum);
				
				else if (coV.covariance)
					// COVARIANCE //
					answer = runCovariance(&results,NULL,NULL,NULL,NULL,NULL,coV.filename1,
										   coV.filename2,one.getMean(),two.getMean(),false,coV,idum);
				else
				{
					cout << "We cannot use covariance and lsf for fitting" << endl;
					cout << "Please specify one or the other. Exiting fit" << endl;
					lateral = false;
				}
				
				cout << setprecision(16) << answer << endl;
				
				if (i == 0) left = answer;
				else if (i == 1) middle = answer;
				else if (i == 2) right = answer;
			}
			
			// Check if results are better than previous
			if ((middle > left) || (middle > right))
			{
				// Check if lower is better
				if (left > right) total_x += one_percent_x;
				else total_x -= one_percent_x;
				
				vertical = true;
			}
			else 
				// We have reached the minimum
			{
				lateral= false;
				cout << "X Translation of Image is at a minimum : " << total_x << endl;
			}
		}
		
		left = right = middle = 0;
		
		if (vertical)
		{
			// VERTICAL SEARCH //
			for (int i = 0; i < 3; i++)
			{
				coV.translate_y = co.translate_y + total_y + ((i-1)*one_percent_y);
				two.readFile(coV.filename2,coV);
				//two.translate(total_x, total_y + ((i-1)*one_percent_y));
				
				// SET RESULTS VARIABLES //
				setComparativeData(counter,coV.filename1,&one,coV.filename2,&two,&results,coV);
				/////////////////////////////////////
				
				if (coV.lsf)
					// LSF // 
					answer = runLsf(&results,NULL,NULL,NULL,NULL,NULL,coV.filename1,coV.filename2,NULL,(int) one.getM(),
									(int) one.getN(),one.getStep(),false,coV,idum);
				
				else if (coV.covariance)
					// COVARIANCE //
					answer = runCovariance(&results,NULL,NULL,NULL,NULL,NULL,coV.filename1,
										   coV.filename2,one.getMean(),two.getMean(),false,coV,idum);
				else
				{
					cout << "We cannot use covariance and lsf for fitting" << endl;
					cout << "Please specify one or the other. Exiting fit" << endl;
					vertical = false;
				}
				
				cout << setprecision(16) << answer << endl;
				
				if (i == 0) left = answer;
				else if (i == 1) middle = answer;
				else if (i == 2) right = answer;
			}
			
			// Check if results are better than previous
			if ((middle > left) || (middle > right))
			{
				// Check if lower is better
				if (left > right) total_y += one_percent_y;
				else total_y -= one_percent_y;
				
				lateral = true;
			}
			else 
				// We have reached the minimum
			{
				vertical = false;
				cout << "Y Translation of Image is at a minimum : " << total_y << endl;
			}
		} 
		
		cout << total_x << " " << total_y << endl;
		
	}
	
	return EXIT_SUCCESS;
}
**/ 
