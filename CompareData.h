#ifndef COMPARE_DATA_H
#define COMPARE_DATA_H

#include <cstdlib>
#include <vector>
#include <math.h>
#include <iostream>

class Compare_Data{
	
public:
	Compare_Data();
	~Compare_Data() {;}
	
	Compare_Data(std::vector<std::string> *o)
	{output_content = o;}
	
	double lsf(const bool d);
	
	double covariance();
	
	void setDataArrayOne(std::vector<float> one, int a, int b) 
	{dataOne = one;  m = a; n = b;}
	
	void setDataArrayTwo(std::vector<float> two) 
	{dataTwo = two;}
	
	void setCount(int a) 
	{count = a;}
	
	std::vector<float> getLSFDifference();
	
private:
	std::vector<std::string> *output_content;
	
	std::vector<float> dataOne;
	std::vector<float> dataTwo;
	std::vector<float> difference;
	
	int m;
	int n;
	int count;
};

#endif
