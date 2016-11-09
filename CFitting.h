/*
 *  fitting.h
 *  zContrast
 *
 *  Created by Andrew Logsdail on 26/01/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

#ifndef FITTING_H
#define FITTING_H
#include "CompareWorkLink.h"
#include "Structures.h"

class Fitting{
	
public:
	Fitting() {;}
	~Fitting() {;}
	
	// Linear fitting methods.
	bool gaussian(int *s, const CreateV &crV, const CompareV &co, std::vector<std::string> *output_content);
	
	bool scale(int *s, const CreateV &crV, const CompareV &co, std::vector<std::string> *output_content);
	
	// Currently not working
	// Two variable fitting methods
	// bool translate(int *s, const CreateV &crV, const CompareV &co);
	
private:
	void printConvergenceError(const int &c, std::vector<std::string> *output_content);

};

#endif

