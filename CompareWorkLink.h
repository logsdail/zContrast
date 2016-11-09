/*
 *  compare_work_link.h
 *  zContrast
 *
 *  Created by Andrew Logsdail on 30/01/2011.
 *  Copyright 2011 University of Birmingham. All rights reserved.
 *
 */

#ifndef __COMPARE_WORK_LINK_HPP__
#define __COMPARE_WORK_LINK_HPP__

#include <cstdlib>
#include "CZContrast.h"
#include "Information.h"
#include "CompareData.h"
#include "Search.h"
#include "CVariables.h"

Information createM(Variables &tempVar, int *seed, std::vector<std::string> *output_content);

Information createM(std::string &image_oufilename, const std::string &structure_filename, const float &g_iX, 
				   const float &g_iY, const float &g_fGridSize, const float &g_fExponent, const float &g_fAlpha, 
				   const float &g_fScaler, const float &g_fNoise, float const &rotateX, float const &rotateY, 
				   float const &rotateZ, bool const &bFourier, bool const &bCrossSection, bool const &bScreen, 
				   bool const &bRotate, const bool &bSave, const std::vector<std::string> &a_elementName, 
					const std::vector<int> &a_atomicNumber, const std::vector<float> &a_atomicRadii, int *seed, std::vector<std::string> *output_content);

void scaleZValues(Information *file1, Information *file2, CompareV coV, std::vector<std::string> *output_content);

double runLsf(Compare_Data *results, std::vector<std::string> *lsfResults, Search *search, const float x,
			  const float y, const float z, const std::string &filename1, const std::string &filename2, const int counter,
			  const int m, const int n, const int step, bool check, CompareV coV, int *seed, std::vector<std::string> *output_content);

double runCovariance(Compare_Data *results, std::vector<std::string> *covarianceResults, Search *search, const float x,
					 const float y, const float z, const std::string &filename1,const std::string &filename2,
					 const float file1_mean, const float file2_mean, bool check, CompareV coV, int *seed, std::vector<std::string> *output_content);

Information checkGridsMatch(bool check, Search *search, Information file1, Information *file2,
							Variables tv, Variables *var, const bool fitting, int *seed, std::vector<std::string> *output_content);

void setComparativeData(const int &counter, const std::string &filename1, Information *file1, const std::string &filename2,
						Information *file2, Compare_Data *results, CompareV coV, std::vector<std::string> *output_content);

float printRegenerateGrid(Information* file1, Information *file2, std::vector<std::string> *output_content);

#endif
