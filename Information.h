#ifndef INFORMATION_H
#define INFORMATION_H

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 **/

// #include <cstdlib>
// #include <vector>
// #include <math.h>
#include "Utils.h"
// #include "Structures.h"
#include <fstream>
#include <iostream>

/**
 @author Andrew Logsdail
 */

// This is going to hold all our information about the dataset we have, and work out averages etc. for us
class Information{
	
	// This is used to store the point values as it is read in
	struct point
	{
		double x;
		double y;
		double z; // Intensity
	};
	
public:
	Information() {idum = NULL; output_content = NULL;}
	~Information() {emptyVectors();}
	
	Information(int *s, std::vector<std::string> *o)
	{ init(s,o); }
	
	void init(int *s, std::vector<std::string> *o);
	
	bool readFile(const std::string &name);
	
	bool readFile(const std::string &file, const CompareV &coV);
	
	bool saveImage(std::string sz);
	
	bool saveImage(const std::string sz, std::vector<float> newArray, const int x, const int y, const float grid);
	
	bool saveImageMatrix(std::string sz);
	
	bool saveImageMatrix(const std::string &sz, std::vector<float> &newArray, const int &x, const int &y, const float &grid);
	
	void reset();
	
	void scale(float scalar);
	
	int getCount() 
	{return count;}
	
	float getXtotal() 
	{return x_total;}
	
	float getYtotal() 
	{return y_total;}
	
	float getZtotal() 
	{return z_total;}
	
	float getMaxX() 
	{return max_x;}
	
	float getMaxY() 
	{return max_y;}
	
	float getMaxZ() 
	{return max_z;}
	
	float getMinX() 
	{return min_x;}
	
	float getMinY() 
	{return min_y;}
	
	float getMinZ() 
	{return min_z;}
	
	float getStep() 
	{return step;}
	
	float getM() 
	{return m;} //
	
	float getN() 
	{return n;} // These are both int but we return floats for calculations
	
	std::string getFilename() 
	{return filename;}
	
	void setFilename(std::string name) 
	{filename = name;}
	
	std::vector<float> getDataArray() 
	{return dataArray;}
	
	// This next one is going to be important for the search routines and comparing each result to the current best
	void setDataArray(std::vector<float> &newArray, const int &x, const int &y, const float &grid) 
	{reset(); setM(x); setN(y); setStep(grid); dataArray = newArray; inputAnalysis();}
	
	void setM(int x_spread)
	{m = x_spread;}
	
	void setN(int y_spread) 
	{n = y_spread;}
	
	void setStep(float gridSize) 
	{step = gridSize;}
	
	float getMean() 
	{return z_total/count;}
	
	// Translate Data
	void translate(const float diff_x, const float diff_y);
	
private:
	CompareV compare_variables;
	
	std::vector<point> data_points, com_x, com_y;
	
	std::vector<float> dataArray;
	
	int count, m, n;
	
	int *idum; // Pointer to random number seed
	std::vector<std::string> *output_content; // Output array
	
	float x_total, y_total, z_total;
	
	float max_x, max_y, max_z, min_x, min_y, min_z;
	
	float step;
	
	std::string filename;
	
	void calcCOM(); // Get COM
	
	bool readZcon(const std::string &name); // Read Zcon
	
	bool readTxt(const std::string &name); // Read Txt

	void inputAnalysis(); // Finds out information when new dataArray is defined. Used for set data
	
	void checkLimits(point tPB);
	
	void calculateDataArray(); // Calculate Array from Data points
	
	void emptyVectors() 
	{data_points.clear(); dataArray.clear(); com_x.clear(); com_y.clear();} // Clean up
	
	void clearCounters();
};

#endif
