#ifndef __CCONTRASTIMAGE__
#define __CCONTRASTIMAGE__

#include <vector>
#include <string>
#include <cstdlib>
#include "Information.h"
#include "fftw3.h"
#include "math.h"

class CContrastImage ///A class which holds the image information essentially a 2D dynamic array
{
public:
	CContrastImage() {;}
	~CContrastImage() {;}
	
	CContrastImage(std::string sz, std::vector<std::string> *o); //Read File
	
	CContrastImage(int nRows, int nColumns, float gridsize, std::vector<std::string> *o); // Set rows, columns, and spacing between points
	
	double getIntensity(int y, int x); // Return the image intensity for a given x y coordinate
	
	void setIntensity(int y, int x, double value);
	
	Information getInformation(const std::string &filename); // Returns data
	
	int getHeight()
	{return m_iRows;} // Return image height
	
	int getWidth()
	{return m_iColumns;} // Return image width
	
	float getGridSize()
	{return m_fGridSize;} // Return grid size
	
	void fourier(float nScale); //Compute fourier transform
	
	// THESE USE aINFORMATION. MUST CLEAR DATA FROM IT WHEN USED //
	void readImage(std::string sz); 
	
	void setDataArray(std::vector<float> dataArray) 
	{m_dArray = dataArray;}
	
  	void setWidth(int x_spread) 
	{m_iColumns = x_spread;}
	
  	void setHeight(int y_spread) 
	{m_iRows = y_spread;}
	
  	void setStep(float gridSize) 
	{m_fGridSize = gridSize;}
	
	void saveImage(std::string sz);
	
	void saveImageMatrix(std::string sz);
	
	void saveCrossSection(const std::string &sz,const int &n);
	
private:
	int m_iRows; 
	
	int m_iColumns;
	
	float m_fGridSize;
	
	std::vector<float> m_dArray;
	
	std::vector<std::string> *output_content;
};

#endif
