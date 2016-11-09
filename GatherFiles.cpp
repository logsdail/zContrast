#include "GatherFiles.h"

/**
 24/08/2011
 - Removed need for bGlobal and bOpt in listFiles()
 - Removed unnecessary padding in lsf and cov output files
 - Allowed the identification of various xyz files, as we can now search through them
 - Removed check on filename for lsf and cov output, as guaranteed to be > 0
 02/10/2011
 - Took out dereferencing of print vector
 24/10/2011
 - Set her up to deal with folders inside folders. Previously this just returned null
 **/

using namespace std;

GatherFiles::GatherFiles() 
{
	// counter = 1; // Appenditure Counter
	lsf_filename = "";
	covariance_filename = ""; // Filenames
}

vector<string> GatherFiles::listFiles( const string &name, const bool &bXYZ)
// This will discover all files within a working directory that we can review
// Inputs: string(dirName) - Directory name
{
	// INPUT VARIABLES //
	const string directory_word="/";
	string sub = name.substr(name.size()-1);
	vector<string> files;
	/////////////////////
	
	// CHECK FOR DIRECTORIES //
	if (cmpStr(sub,directory_word))
	{
		if (bXYZ) 
		{
			files = listXYZFiles(name.substr(0,name.size()-1));
		}
		else 
		{
			files = listImageFiles(name.substr(0,name.size()-1));
		}
	}		
	else 
		files.push_back(name);
	
	return files;
}

vector<string> GatherFiles::listImageFiles( const string &dirName )
// This will discover all files within a working directory that we can review
// Inputs: string(dirName) - Directory name
{
	vector<string> files;
	
	DIR *dirp = opendir( dirName.c_str() );
	if (dirp)
	{
		struct dirent *dp = NULL;
		while ((dp = readdir( dirp )) != NULL )
		{
			string file( dp->d_name );
			
			// cout << file << endl;
			
			if ( file == "." || file == ".." )
			{
				// skip these
				continue;
			}
			
			if ( dp->d_type & DT_DIR )
			{
				// found a directory; recurse into it.
				string filePath = dirName + "/" + file;
				vector<string> files_temp = listImageFiles(filePath);
				
				files.insert(files.end(), files_temp.begin(), files_temp.end());
			}
			else
			{
				//regular file found				
				const string zcon="zcon";
				const string txt=".txt";
				string sub = file.substr(file.size()-4);
				
				// cout << file << endl;
				
				if (cmpStr(sub,zcon) || cmpStr(sub,txt))
				{
					// cout << sub << endl;
					files.push_back(dirName + "/" + file);
				}
			}
		}
		closedir(dirp);
	}
	
	return files;
}

vector<string> GatherFiles::listXYZFiles( const string &dirName )
// This will discover all files within a working directory that we can review
// Inputs: string(dirName) - Directory name
{
	vector<string> files;
	
	DIR *dirp = opendir( dirName.c_str() );
	if (dirp)
	{
		struct dirent *dp = NULL;
		while ((dp = readdir( dirp )) != NULL )
		{
			string file( dp->d_name );
			
			if ( file == "." || file == ".." ) 
			{	
				// skip these
				continue;
			}
			
			if ( dp->d_type & DT_DIR )
			{
				// found a directory; recurse into it.
				string filePath = dirName + "/" + file;
				vector<string> files_temp = listXYZFiles(filePath);
				files.insert(files.end(), files_temp.begin(), files_temp.end());
			}
			else
			{
				//regular file found
				const string xyz=".xyz";
				string sub = file.substr(file.size()-4);
				
				if (cmpStr(sub,xyz))
				{
					files.push_back(dirName + "/" + file);
				}
			}
		}
		closedir(dirp);
	}
	
	return files;
}

bool GatherFiles::prepareOutputFiles(bool check, bool bLsf, const string lsf, bool bCov, const string cov)
// Prepare the output files for data to be written to them
// Inputs - check to see if we are saving data
//		  - check if we are running lsf analysis
//		  - lsf filename
//		  - check if we are running covariance analysis
//		  - covariance filename
{
	// Error check
	bool error = false;
	
	if (check)
	{
		vector<string> empty;
		
		//cout << bLsf << " " << lsf << endl;
		//cout << bCov << " " << cov << endl;
		
		if (bLsf && !error) 
		{
			setLsfFilename(lsf);
			error = write(empty,lsf);
		}
		
		if (bCov && !error)
		{
			setCovarianceFilename(cov);
			error = write(empty,cov);
			// cout << cov << " " << covariance_filename << endl;
		}
	}
	
	// cout << lsf_filename << " " << covariance_filename << endl;
	
	// Return error check
	return error;
}

bool GatherFiles::write(vector<string> data, const string filename)
// This is just a text outputting function to file
// Inputs vector<string>(data) - Data string to be written
// 	  string(filename) - File to be written too
// Outputs: Success or Failure
{
	ofstream outData;
	outData.open(filename.c_str(), ios::out);
  	if(!outData.is_open())
	{
		cout << "Error writing to file: " << filename << endl;
		return true;
	}
	else
	{
		//outData << "=== Data Output to : " << filename << " === " << endl;
		
		for (vector<string>::iterator it = data.begin(); it != data.end(); ++it)
		{
			outData << *it << endl;
		}
		outData.close();
		return false;
	}
}

bool GatherFiles::appendToOutputFiles(vector<string> lsf, vector<string> cov)
// Append data to output files
// Inputs - lsf data
//		  - covariance data
{
	// Error check
	bool error = false;
	// cout << lsf_filename.length() << " " << lsf_filename << " " << covariance_filename.length() << " " << covariance_filename << endl;
	
	if (lsf_filename.size() > 0)
	{
		error = append(lsf,lsf_filename);
	}
	
	if (covariance_filename.size() > 0)
	{
		append(cov,covariance_filename);
	}
	
	// Return error check
	return error;
}

bool GatherFiles::append(vector<string> data, const string filename)
// This is just a text appending function to file
// Inputs vector<string>(data) - Data string to be written
//        string(filename) - File to be written too
// Outputs: Success or Failure
{
	ofstream outData;
	outData.open(filename.c_str(), ios::app);
	
	if(!outData.is_open())
	{
		cout << "Error appending to file: " << filename << endl;
		return false;
	}
	else
	{
	
		
		for (vector<string>::iterator it = data.begin(); it != data.end(); ++it)
		{
			outData << *it << endl;
		}
		outData.close();
		
		//counter++;
		return true;
	}
}
