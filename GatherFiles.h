#ifndef GATHERFILES_H
#define GATHERFILES_H

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 **/

// #include <cstdlib>
// #include <vector>
#include <iostream>
#include <fstream>
#include "dirent.h"
#include "sys/types.h"
#include "Utils.h"

/**
 @author Andrew Logsdail
 */

// I am using this class to list all files in a directory
// as well as write text outputs
// Currently data outputs are written through aInformation.cpp

/**
 24/08/2011
 - Updated listFiles to reflect changes in cpp
 - Removed counter for appendages
 **/

class GatherFiles{
	
public:
	GatherFiles();
	~GatherFiles() {;}
	
	//std::vector<std::string> listFiles(const std::string &name, const bool &bOpt, const bool &bGlobal, const bool &bBoth);
	std::vector<std::string> listFiles(const std::string &name, const bool &bBoth);
	
	bool prepareOutputFiles(bool check, bool bLsf, const std::string lsf, bool bCov, const std::string cov);
	
	bool write(std::vector<std::string> data, const std::string filename);
	
	bool appendToOutputFiles(std::vector<std::string> lsf, std::vector<std::string> cov);
	
	bool append(std::vector<std::string> data, const std::string filename);
	
	void setLsfFilename(std::string f) 
	{lsf_filename = f;}
	
	void setCovarianceFilename(std::string f) 
	{covariance_filename = f;}
	
private:
	// int counter;
	
	std::string lsf_filename;
	std::string covariance_filename;
	
	std::vector<std::string> listXYZFiles(const std::string &dirName);
	
	std::vector<std::string> listImageFiles(const std::string &dirName);
	
};

#endif
