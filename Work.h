#ifndef WORK_H
#define WORK_H

#include <vector>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include "CompareWorkLink.h" 
#include "GatherFiles.h"
#include "ReadVariables.h"
#include "Structures.h"
#include "CFitting.h"

class Work{
	
public:
	Work();
	~Work() {;}
	
	void setVariables(ReadVariables re);
	
	void setProcessors(const int p);
	
	void setSeed(int s) 
	{ seed = s; }
	
	void setJob(std::string job) 
	{ function = job; }
	
	void run();
	
	void create();
	
	void fit();
	
	//bool search();
	
	//bool search_MPI(const int numprocs);
	void search_MPI();
	
	void slave_MPI(const std::string filename1, const std::string filename2, const int myrank,
				   const bool bBoth, const bool bSearch, const bool bOpt, const bool bGlobal);
	
private:
	int seed;
	int myrank;
	// Setup output files
	std::string output_log;					 // Output filename //
	std::string tmpdir;
	std::string pwd;
	
	std::vector<std::string> output_content; // Output log //
	std::string function;
	
	ReadVariables red;
	
	//std::vector<Information> preloadImages(std::vector<std::string> f);
	
	void printErrorDifferentSizes(const std::string &filename1, const int &f1count, const std::string &filename2, const int &f2count);
	
	void printFileReadError(const std::string &fn);
	
	void printFileLoaded(const std::string &fn);
	
	void printComparisons(const std::string &filename1, const std::string &filename2, const int &c, const int &cc, RotationPoint min);
};
#endif
