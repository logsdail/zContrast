#include "Work.h"

/**
 24/08/2011
 - Moved second file loop outside parallelisation. This means we can run for more than one type of xyz file
 a GA search, but that cannot parallelise straight comparisons. Not seen as a drawback as straight comparisons
 usually don't need OMP anyway
 - Redefined output files for lsf and cov, so now they reflect filename
 - Added functions to restart search, as need for multiple comparisons.
 - Moved some variable declarations outside for loops.
 - Print comparison writes file names, otherwise we won't know on parallelised outputs.
 - Need to change to MPI. Good examples of what to do at:
 http://sgowtham.net/blog/2010/11/28/mpi-c-an-advanced-send-receive/
 http://www.lam-mpi.org/tutorials/one-step/ezstart.php
 Data structures:
 http://www.lam-mpi.org/tutorials/one-step/datatypes.php
 https://computing.llnl.gov/tutorials/mpi/#Derived_Data_Types
 Hybrid OMP/MPI:
 http://www.slac.stanford.edu/comp/unix/farm/mpi_and_openmp.html
 25/08/2011
 - Removed image file preloading. This is not time consuming and now we can parallelise easier.
 - Incorporated MPI. Set up Search_MPI() and Slave_MPI() methods. Need to consolidate the content
 26/08/2011
 - Consolidated Search method.
 - Implemented printscreen output through a string vector, written to file every full complete of search
 14/09/2011
 - Set up outputs to write to temp directory, then copy the data back to working directory at the end.
 This saves on cross network IOs during the run. Limitation is we can't use more than one node for each GA search,
 but this isn't currently implementated anyway so is not a big problem.
 30/09/2011
 - Fixed bugs to do with unassigned pointers giving malloc errors. Need to reimplement the writing of IOs to temp directory,
 as this was disabled as part of the debugging.
 - Corrected so now this works fine - it writes to tmp in all occassions.
 - In the constructor we also check if the directory labels end in slashes, to avoid miswritten files.
 02/10/2011
 - Set up comparative data files to also be copied backwards and forwards from tmp
 24/10/2011
 - Pumps out myrank to check MPI config
 28/10/2011
 - Updated MPI to workload using the master node as well as slaves
 - Edited readTxt to work successfully so we can load a txt file and output the zcon of it.
 **/

using namespace std;

Work::Work() 
// Constructor - set seed to 0
{
	seed = 0;
	//** INTIALISE MPI **//
	myrank = 0;
	
	// Setup output files
	NumberToString(myrank,output_log);
	output_log = "process" + output_log + ".log";
	
	// Need to check if these end in slashes or not
	tmpdir = getenv("TMPDIR");
	unsigned int found = tmpdir.find_last_of('/');
	if (found != (tmpdir.length() - 1))
	{
		tmpdir += "/";
	}
	
	pwd = getenv("PWD");
	found = pwd.find_last_of('/');
	if (found != (pwd.length() - 1))
	{
		pwd += "/";
	}
}

void Work::setProcessors(const int p)
// Set number of processors for parallel runs
// Inputs(int) - number of processors
{
	omp_set_num_threads(p);
}

void Work::setVariables(ReadVariables re)
// Store all variables used in the calculation here so we can use them as needed
// Input(re) - All variables for program
{
	red = re;
	
	// Set up all variables specific to this system
	setProcessors(red.ma_v.processors);
	setSeed(-red.ma_v.seed);
	setJob(red.ma_v.function);
}

void Work::run()
// Designed to run the program. Checks to see what we are doing.
{		
	GatherFiles directory;

	if (cmpStr(red.ma_v.function,"create"))
	{
		// Not parallelised so we'll just write to working directory
		create();
		directory.write(output_content,pwd+output_log);
	}
	else if (cmpStr(red.ma_v.function.substr(0,4),"fit_"))
	{
		// Not parallelised so we'll just write to working directory
		fit();
		directory.write(output_content,pwd+output_log);
	}
	else 
	{
		search_MPI();
	}
}

void Work::create()
//bool Work::create(const string &filename)
// Creates new Zcontrast images using data from input file
{
	Variables var(red.getCreateVariables().xyz_filename,&seed);
	var.setVariables(red.getCreateVariables());
	var.bSave = true;

	// START STRUCTURE CALCULATION //
	Information check = createM(var,&seed,&output_content);
	/////////////////////////////////
} 

void Work::fit()
{
	Fitting fit;
	//fit.setVariables(&seed,red.getCreateVariables(),red.getCompareVariables());
	
	if (cmpStr(function.substr(4,function.length()-4),"gaussian"))
	{
		fit.gaussian(&seed,red.getCreateVariables(),red.getCompareVariables(),&output_content);
	}
	else if (cmpStr(function.substr(4,function.length()-4),"scale"))
	{
		fit.scale(&seed,red.getCreateVariables(),red.getCompareVariables(),&output_content);
	}
	//else if (cmpStr(function.substr(4,function.length()-4),"translate"))
	//{
	//	fit.translate(&seed,red.getCreateVariables(),red.getCompareVariables());
	//}
	else 
	{
		cout << "No fitting function for >> " << function.substr(4,function.length()-4) << endl;
	}
	
	
}

void Work::search_MPI()
//bool Work::search_MPI(const int numprocs)
// Loads data, creates images where necessary and compares it using compare class functions
{
	int complete = 0;
	int t_complete = 0;
	int success = 1;
	unsigned int terminate = LARGE;
	const int master = 0;
	MPI_Status status;
	// Check if we are running MPI
	int numprocs = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	if (numprocs > 1)
	{
		int namelen;
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		
		MPI_Get_processor_name(processor_name, &namelen);
		// MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		
		// Get TMPDIR system variable
		// We will then write to this to save cross-network IOs
		// And copy the data back to the user area at the end
		// string tmpdir = getenv("TMPDIR");
		// Setup output files
		NumberToString(myrank,output_log);
		output_log = "process" + output_log + ".log";
		setSeed(seed-myrank); // Sets individual random number seed for each process
		
		if (myrank == 0) cout << "Master speaking: We are running in MPI mode" << endl;
		cout << "This is process: " << myrank << " out of: " << numprocs << ". Processor name: " << processor_name << endl;
	}
	else
	{
		cout << "We are not running in MPI mode. Numprocs = " << numprocs << ", Myrank = " << myrank << endl;
	}
	
	// CHECK IF WE ARE JUST COMPARING OR DOING BOTH CREATE AND COMPARE //
	// This is necessary for input file analysis
	bool bBoth = false;
	bool bSearch = false;
	bool bOpt = false; // OPTIMISATION
	bool bGlobal = false; // Global Process
	
	if (cmpStr(red.ma_v.function,"both")) 
	{
		bBoth = true;
	}
	else if (cmpStr(red.ma_v.function,"search")) 
	{
		bBoth = true;
		bSearch = true;
		
		bOpt = red.getSearchVariables().bMinimise;
		bGlobal = red.getSearchVariables().bGlobal_search;
	}
	
	// INITIALISE File VARIABLES //
	GatherFiles directory;
	string filename1 = red.getCompareVariables().filename1;
	string filename2 = red.getCompareVariables().filename2;
	vector<string> fileList1, fileList2;
	fileList1 = directory.listFiles(filename1,bBoth);
	fileList2 = directory.listFiles(filename2,false);
	
	vector<string> prep;
	// prep.push_back("");
	directory.write(prep,tmpdir+output_log);
	
	cout << "Files in folder 1: " << fileList1.size();  
	cout << ". Files in folder 2: " << fileList2.size() << endl;
	
	unsigned int count1 = 0;
	unsigned int count2 = 0;
	int max_count = fileList1.size() * fileList2.size();
	
	if (max_count > 1)
	{
		if (numprocs > 1)
		{
			if (myrank == 0)
				// I am the master
				// Let's figure out the filenames and distribute them.
			{
				if (numprocs > max_count)
				{
					// Quick error message in case we have given ourselves too many MPI processors
					cout << "We don't need this many processors in the communicator! Currently: " << numprocs << endl;
					cout << "As we only have " << max_count << " calculations to perform" << endl;
				}
				
				for (int rank = 1; rank < numprocs; ++rank)
				{
					MPI_Send(&count1,1,MPI_INT,rank,0,MPI_COMM_WORLD);
					MPI_Send(&count2,1,MPI_INT,rank,0,MPI_COMM_WORLD);

					count2++;
					
					if (count2 >= fileList2.size())
					{
						count2 = 0;
						count1++;
					}
					
					if (count1 >= fileList1.size())
					{
						// Escape route from too many processors
						rank = numprocs;
					}
				}
				
				cout << "Processes all have jobs" << endl;
				
				// If we have outstanding jobs keep collecting them and redistributing
				while (count1 < fileList1.size())
				{
					// Check how many jobs we have sent out. If everything has a job run the next one on the master
					if (((count1*fileList2.size())+count2+1)%numprocs == 0)
					{
                                                // Now we need to work....
                                                slave_MPI(fileList1[count1],fileList2[count2],myrank,bBoth,bSearch,bOpt,bGlobal);
						t_complete = 1;
					}
					// Otherwise distribute out some more jobs
                                        else
					{
	                                        MPI_Recv(&t_complete,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);

                                        	MPI_Send(&count1,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
                                        	MPI_Send(&count2,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
					}

					//MPI_Recv(&t_complete,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
					complete += t_complete; // Keep count of number left to compute
					
					//MPI_Send(&count1,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
					//MPI_Send(&count2,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
					
					count2++;
					
					if (count2 >= fileList2.size())
					{
						count2 = 0;
						count1++;
					}
				}
				
				// Receive all outstanding jobs
				while (complete < max_count)
				{
					MPI_Recv(&t_complete,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
					complete += t_complete;
				}
				
				cout << "Sending termination signal" << endl;
				// Send out termination signal
				for (int rank = 1; rank < numprocs; ++rank)
				{
					MPI_Send(&terminate,1,MPI_INT,rank,0,MPI_COMM_WORLD);
					MPI_Send(&terminate,1,MPI_INT,rank,0,MPI_COMM_WORLD);
				}
			}
			else
			{
				// Here we will receive the filenames
				MPI_Recv(&count1,1,MPI_INT,master,0,MPI_COMM_WORLD,&status);
				MPI_Recv(&count2,1,MPI_INT,master,0,MPI_COMM_WORLD,&status);
				
				while ((count1 != terminate) && (count2 != terminate))
				{
					// Now we need to work....
					slave_MPI(fileList1[count1],fileList2[count2],myrank,bBoth,bSearch,bOpt,bGlobal);
					
					// And then return confirmation we have finished
					MPI_Send(&success,1,MPI_INT,master,0,MPI_COMM_WORLD);
					
					// And receive new filenames
					MPI_Recv(&count1,1,MPI_INT,master,0,MPI_COMM_WORLD,&status);
					MPI_Recv(&count2,1,MPI_INT,master,0,MPI_COMM_WORLD,&status);
				}
			}
			
			// Get CURRENT DIRECTORY system variable
			// We will then write use this to save cross-network IOs
			// And copy the data back to the user area
			// string pwd = getenv("PWD");
			// string process_log = "";
			// Copy files back from tmps
			// NumberToString(myrank,process_log);
			// process_log = pwd + "/process" + process_log + ".log";
			
			// Run system command to copy files back
			///if (myrank == 0)
			// {
			// 	cout << "Copying log files back to working directory" << endl;
			// }
			// string cmd = "cp " + output_log + " " + process_log;
			// system(cmd.c_str());
		}
		else
			// Not running in MPI mode so just do a normal run
		{
			for (count1 = 0; count1 < fileList1.size(); count1++)
			{
				for (count2 = 0; count2 < fileList2.size(); count2++)
				{
					slave_MPI(fileList1[count1],fileList2[count2],myrank,bBoth,bSearch,bOpt,bGlobal);
				}
			}
		}
	}
	// Otherwise we only have one case to review so no need for MPI
	else
		slave_MPI(fileList1[count1],fileList2[count2],myrank,bBoth,bSearch,bOpt,bGlobal);
	
	// Run system command to copy files back
	string cmd = "cp " + tmpdir + output_log + " " + pwd + output_log;
	// cout << "Copying log files back to working directory : " << cmd << endl;
	system(cmd.c_str());
	
	// return EXIT_SUCCESS;
}

void Work::slave_MPI(const string filename1, const string filename2, const int myrank,
					 const bool bBoth, const bool bSearch, const bool bOpt, const bool bGlobal)
{
	// Print to screen output. Now we will try to pipe all print screens to a string, and append them to output_log
	cout << "Running on processor " << myrank << " (Log file is " << tmpdir << output_log << ") with files: " << filename1 << " and " << filename2 << endl;
	
#pragma omp critical
	{
		output_content.push_back("=================================================");
		output_content.push_back("Searching with the files: " + filename1 + " and " + filename2);
		output_content.push_back("=================================================");
	}
	///////////
	Search search;
	Variables var(&seed);
	
	if (bBoth)
	{
		var.setVariables(red.getCreateVariables());
	}
	
	// RESULTS VARIABLES //
	Compare_Data results(&output_content);
	vector<string> lsfResults; // LSF //
	vector<string> covarianceResults; // VARIANCE //
	string output_lsf;
	string output_cov;
	string f1;
	/////////////////////////
	
	// DATA VARIABLES//
	Information file1(&seed,&output_content);
	Information file2(&seed,&output_content);
	file2.readFile(filename2,red.getCompareVariables());
	
	// INTEGER COUNTERS
	int counter = 0;
	int pos = 0;
	int i = 0;
	int loops_this_time = 0; 
	
	// CLEAR RESULTS FILES //
	pos = filename1.find_last_of('/');
	// output_lsf = f1;
	output_lsf = filename1.substr(pos+1) + "_";
	// output_cov = f1;
	output_cov = filename1.substr(pos+1) + "_";
	
	pos = filename2.find_last_of('/');
	output_lsf += filename2.substr(pos+1) + ".lsf";
	output_cov += filename2.substr(pos+1) + ".cov";
	
	GatherFiles directory;
	directory.prepareOutputFiles(red.getCompareVariables().save_results,red.getCompareVariables().lsf,
								 tmpdir+output_lsf, red.getCompareVariables().covariance, tmpdir+output_cov);
	
	// SETUP FOR SEARCH //
	search.init(&seed,&output_content);
	if (bSearch)
	{
		search.getVariables(red.getSearchVariables(),var);
	}
	//////////////////////
	
	// WHILE LOOP TO SEE IF WE HAVE FINISHED SEARCHING //
	while (!search.getComplete())
	{
		// IF THIS IS A SEARCH SET LOOPS FOR FOR LOOP, ELSE WE'LL JUST GO ONCE //
		if (bSearch) 
			loops_this_time = search.getRotationsArraySize();
		else
		{
			search.setComplete();
			loops_this_time = 1;
		}
		
		////////////////////////////////////////////////////////////////////////
		// PARALLELISE FROM HERE //
		///////////////////////////
		
		// Original working OMP constructs on first line
#pragma omp parallel for default(none) \
private(f1,file1,results,i) \
shared(file2,filename1,filename2) \
shared(loops_this_time,counter) \
shared(search,var) \
shared(lsfResults,covarianceResults,cout) 
		// Original working OMP constructs
		
		for (i = 0; i < loops_this_time; i++)
		{
			// cout << i << " : " << loops_this_time << endl;
			
			Variables tv = var;
			f1 = filename1;
			
			// GET CREATE DATA FOR LATER COMPARISON, ELSE JUST READ FILE //
			if (bBoth)
			{
				//IF WE ARE SEARCHING GET NEXT ATTEMPT AT SEARCH //
				// Critical should remove any remaining seg. faults
				if (bSearch) 
#pragma omp critical
				{
					tv = search.getNext(i);
				}
				
				//SET FILENAME AND POINTERS //
				tv.init(f1,&seed);
				// tv.setStructureFilename(f1);
				//CREATE STRUCTION //
				file1 = createM(tv,&seed,&output_content);
				//GET NEW FILENAME WITH ROTATIONS //
				f1 = file1.getFilename();
			}
			else 
			{
				file1.readFile(f1,red.getCompareVariables());
			}
			//////////////////////////////////////////////////////////////
			
			if (file1.getCount() != 0)
			{
				printFileLoaded(f1);
				
				if (file2.getCount() != 0)
				{
					printFileLoaded(filename2);
					
					// Check images are the same size
					file1 = checkGridsMatch((bBoth && (file1.getCount() != file2.getCount())), &search, file1,
											&file2, tv, &var, false, &seed, &output_content);
					
					// CHECK IF WE HAVE FILES OF MATCHING SIZE //
					if ((file1.getCount() == file2.getCount()))
					{
						// If we have got this far we have two sets of data, of equal size
						counter++; //In case we need to know how many comparisons we have done
						
						// SCALE RESULTS //
						scaleZValues(&file1,&file2,red.getCompareVariables(),&output_content);
						
						// SET RESULTS VARIABLES //
						setComparativeData(counter,f1,&file1,filename2,&file2,&results,red.getCompareVariables(),&output_content);
						/////////////////////////////////////
						
						double a;
						// LSF //  
						a = runLsf(&results,&lsfResults,&search,tv.rotateX,tv.rotateY,tv.rotateZ,f1,filename2,
								   counter,(int) file1.getM(),(int) file1.getN(),(int) file1.getStep(),
								   ((bOpt || bGlobal) && cmpStr(search.getOptimisationRef(),"lsf")),red.getCompareVariables(),
								   &seed,&output_content);
						
						// COVARIANCE //
						a = runCovariance(&results,&covarianceResults,&search,tv.rotateX,tv.rotateY,tv.rotateZ,f1,
										  filename2,file1.getMean(),file2.getMean(),
										  ((bOpt || bGlobal) && cmpStr(search.getOptimisationRef(),"covariance")),
										  red.getCompareVariables(),&seed,&output_content);
						
						if (red.getCompareVariables().print_to_screen) 
#pragma omp critical
						{
							output_content.push_back("");
						}
					}
					else 
					{
						printErrorDifferentSizes(f1,file1.getCount(),filename2,file2.getCount());
					}
				}
				else 
				{
					printFileReadError(filename2);
				}
			} 
			else 
			{
				printFileReadError(f1);
			}
			
		}
		
		cout << "Updating output file: " << tmpdir+output_log << " with " << output_content.size() << " lines." <<  endl;
		// APPEND RESULTS TO CURRENT DATA //
                cout << "Updating results files. LSF: " << lsfResults.size() << ", Covariance : " << covarianceResults.size() << endl;
		directory.appendToOutputFiles(lsfResults,covarianceResults);
                cout << "Updating output file: " << tmpdir+output_log << " with " << output_content.size() << " lines." <<  endl;
		directory.append(output_content,tmpdir+output_log);
		/////////// CLEAR VECTORS //////////
		lsfResults.clear();        
		covarianceResults.clear();
		output_content.clear();
		////////////////////////////////////
	}
	
	// Copy comparative output files back across
	if (red.getCompareVariables().lsf)
	{
		string cmd = "cp " + tmpdir + output_lsf + " " + pwd + output_lsf;
		system(cmd.c_str());
	}
	
	if (red.getCompareVariables().covariance)
	{
		string cmd = "cp " + tmpdir + output_cov + " " + pwd + output_cov;
		system(cmd.c_str());
	}
	
	// Print Comparisons
	if (bSearch)
	{
		printComparisons(filename1,filename2,counter,search.getHistoryCopyCounter(),search.getMinimumPoint());
	}
}

void Work::printErrorDifferentSizes(const string &filename1, const int &f1count, const string &filename2, const int &f2count)
// Error message different array sizes
{
	string f1countS;
	string f2countS;
	NumberToString(f1count,f1countS);
	NumberToString(f2count,f2countS);
	
#pragma omp critical
	{
		output_content.push_back("Data sets are not of matching size. " + filename1 + ": " + f1countS + ", " + filename2 + ": " + f2countS);
	}
}

void Work::printFileReadError(const string &fn)
// Output file read error
{
#pragma omp critical
	{
		output_content.push_back("Error reading file: " + fn);
	}
}

void Work::printComparisons(const string &filename1, const string &filename2, const int &c, const int &cc, RotationPoint min)
// Print total number of comparisons
{
//#pragma omp critical
//	{
		cout << endl;
		cout << "For " << filename1 << " and " << filename2 << " ... " << endl;
		cout << "A total of " << c << " comparisons (FEs) were made, with a further " << cc << " copied from history." << endl;
		cout << "The found optimised orientation is (" << min.theta << " , " << min.phi << " , " << min.psi << " ) : " << min.value << endl;
		cout << endl;
//	}
	
	string cS;
	string ccS;
	NumberToString(c,cS);
	NumberToString(cc,ccS);
	
#pragma omp critical
	{
		output_content.push_back("For " + filename1 + " and " + filename2 + " ... ");
		output_content.push_back("A total of " + cS + " comparisons (FEs) were made, with a further " + ccS + " copied from history.");
	}
}

void Work::printFileLoaded(const string &fn)
// Print to screen if file is loaded successfully
{
	// SET FILE VALUES AFTER SUCCESFUL LOADING //
	if (red.getCompareVariables().print_to_screen) 
	{
#pragma omp critical
		{
			output_content.push_back("File Loaded : " + fn);
		}
	}
}
