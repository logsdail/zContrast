#include "Information.h"

/**
 24/10/2011
 - Edited COM to move pixels if they are greater than 50% of a step from desired search area - dealt with 
 problems assoicated with moving pixels into random areas of search area.
 - Edited read txt to no longer be reliant on the pointer to a random number seed
 **/

using namespace std;

void Information::init(int *s, vector<string> *o) 
// Constructor - clears all variables ready for use
// Takes seed as input variable
{
	this->reset();
	// Initiate Random Numbers
	idum = s;
	output_content = o;
}	

void Information::reset()
// Used to reset all data to 0.
// Called when initialising data sets, takes no inputs
{
	this->clearCounters();
	this->emptyVectors();
}

void Information::clearCounters()
// Used to reset all Counters, not vectors
{
	x_total = 0.0;
	y_total = 0.0;
	z_total = 0.0;
	
	max_x = 0;
	max_y = 0;
	max_z = 0;
	
	min_x = LARGE;
	min_y = LARGE;
	min_z = LARGE;
	
	count = 0;
	m = 0;
	n = 0;
	step = 0;
}

bool Information::readFile(const string &name)
// Reads in new data
// Input(name) - filename
// Outputs - boolean of success or failure
{
	CompareV c;
	c.filename1 = c.filename2 = "";
	//	c.lsf_file = c.covariance_file = "";
	c.scale_intensity = false;
	c.print_to_screen = false;      
	c.lsf = false;                  
	c.covariance = false;            
	c.save_results = false;        
	c.centre_image = false;
	c.translate = false;
	c.save_lsf_difference = false;
	c.translate_x = 0;
	c.translate_y = 0;
	
	return readFile(name,c);
}

bool Information::readFile(const string &file, const CompareV &coV)
// Reads in new data
// Input(name) - filename 1 or 2
//      (custom) - list of all used boolean operators
// Outputs - boolean of success or failure
{
	reset();
	compare_variables = coV;
	setFilename(file);
	
	string sub = filename.substr(filename.size()-4);
	
	if (cmpStr(sub,"zcon"))
	{
		return readZcon(filename);
	}
	else if (cmpStr(sub,".txt"))
	{
		return readTxt(filename);
	}
	else 
	{
		return false;
	}
}

bool Information::saveImage(string sz)
// Saves data currently held in zcon format
// Input(sz) - filename
// Outputs - boolean of success or failure
{
	ofstream fsFile;
	fsFile.open(sz.c_str(), ios::out);
	int position;
	
	if(!fsFile.is_open())
	{
#pragma omp critical
		{
			output_content->push_back("Error Writing output file: " + sz);
		}
		return EXIT_FAILURE;
	}
	else
	{
		for(int y = 0; y < getN(); y++)
		{
			for(int x = 0; x < getM(); x++)
			{
				position = (y*n) + x;
				fsFile << x * getStep() << "\t" << y * getStep() << "\t" << dataArray[position] << endl;
			}
			fsFile << endl; //we require this for pm3d plot to work in gnuplot
		}
		fsFile.close();
		return EXIT_SUCCESS;
	}
}

bool Information::saveImage(const string sz, vector<float> newArray, const int x, const int y, const float grid)
// Saves data passed to class as Zcon
// Input string(sz) - filename
//	 vector<float>(newArray) - Data to be saved
//	 int(x) - width
//	 int(y) - height
//	 float(grid) - step sizes
// Outputs - boolean success or failure
{
	this->reset();
	dataArray = newArray;
	this->setM(x);
	this->setN(y);
	this->setStep(grid);
	return saveImage(sz);
}

bool Information::saveImageMatrix(string sz)
// Saves data passed to class as Txt
// Input string(sz) - filename
// Outputs - boolean success or failure
{
	ofstream fsFile;
	fsFile.open(sz.c_str(), ios::out);
	int position;
	
	if(!fsFile.is_open())
	{
#pragma omp critical
		{
			output_content->push_back("Error Writing output file: " + sz);
		}
		return EXIT_FAILURE;
	}
	else
	{
		for(signed int y = (int) getN(); y > -1; y--)
		{
			for(int x = 0; x < getM(); x++)
			{
				position = (y*n) + x;
				// cout << position << " " << dataArray.size() << endl;
				fsFile << dataArray[position] << " ";
			}	
			fsFile << std::endl; //we require this for pm3d plot to work in gnuplot
		}
		fsFile.close();
		return EXIT_SUCCESS;
	}
}

bool Information::saveImageMatrix(const string &sz, vector<float> &newArray, const int &x, const int &y, const float &grid)
// Saves data passed to class as Txt matrix
// Input string(sz) - filename
//       vector<float>(newArray) - Data to be saved
//       int(x) - width
//       int(y) - height
//       float(grid) - step sizes
// Outputs - boolean success or failure
{
	this -> reset();
	dataArray = newArray;
	this->setM(x);
	this->setN(y);
	this->setStep(grid);
	return saveImageMatrix(sz);
}


bool Information::readTxt(const string &name)
// Read in text file. Perhaps we should change this to you use the tokenizer at some point
// Inputs: string(name) - filename
// Outputs - boolean success or failure
{
	ifstream inData;
	inData.open( name.c_str() );
	string buffer;
	
	if ( !inData ) 
	{
		return EXIT_FAILURE;
	}
	else
	{
		int lines = 0;
		while ( getline(inData,buffer) )
		{
			if (buffer.size() != 0)
			{
				lines++;
			}
		}
		inData.clear();
		inData.seekg(ios::beg);
		lines--; // We seem to over read by one line
		
		unsigned int i, j, x, y, string_size;
		string ts;
		point tPB,oldPB;
		float dx, odx, dy, ody;
		dx = odx = dy = ody = 0;
		
		//oldPB.x = randomNumber(101,idum);
		//oldPB.y = randomNumber(101,idum);
		//oldPB.z = randomNumber(101,idum);
		oldPB.x = 77777;
		oldPB.y = 88888;
		oldPB.z = 99999;
		y = lines; 
		
		while ( getline(inData,buffer) )
		{
			j = 0;
			i = 0;
			x = 0;
			while( i < buffer.size())
			{
				if(buffer.size() != 0);
				{
					j=i;
					find_character_end(buffer,j);
					string_size = j-i;
					ts = buffer.substr(i,string_size);
					i=j;
					find_white_space_end(buffer,j);
					i=j;
					
					tPB.x = x;
					tPB.y = y;
					StringToNumber(ts,tPB.z);
					
					data_points.push_back(tPB);
					//inputAnalysis(tPB);
					
					odx = dx;
					dx = (x - oldPB.x);
					if ((odx == dx) && (ody == dy))
					{
						step = dx;
					}
					
					oldPB = tPB;
					x++;
				}
			}
			ody = dy;
			dy = (y - oldPB.y);
			y--;
		}	
	}
	inData.close();
	
	inputAnalysis();
	
	if (compare_variables.translate) 
	{	
		if (compare_variables.centre_image) 
		{
			calcCOM(); // Calculate Centre of Mass
		}
		else 
		{
			translate(compare_variables.translate_x,compare_variables.translate_y);
		}
	}
	calculateDataArray();
	
	// Save Centred Image if need be
	if (compare_variables.translate && compare_variables.save_results)
	{
		if (compare_variables.centre_image) 
		{
			saveImage("COM_" + name.substr(0,name.length()-4) + ".zcon");
		}
		else 
		{
			saveImage("TRANS_" + name.substr(0,name.length()-4) + ".zcon");
		}
	}
	
	return EXIT_SUCCESS;
}

bool Information::readZcon(const string &name)
// Read in zcon file. Perhaps we should change this to you use the tokenizer at some point
// Inputs: string(name) - filename
// Outputs - boolean success or failure
{
   	ifstream inData;
   	inData.open( name.c_str() );
	
   	if ( !inData ) 
	{
		return EXIT_FAILURE;
	}
	else 
	{
		point tPB, oldPB;
		tPB.x = oldPB.x = randomNumber(101,idum);
		tPB.y = oldPB.y = randomNumber(101,idum);
		tPB.z = oldPB.z = randomNumber(101,idum);
		float dx, odx, dy, ody;
		dx = odx = dy = ody = 0;
		
   		while (!inData.eof())
   		{
			inData >> tPB.x;
			inData >> tPB.y;              
			inData >> tPB.z;
			
			// cout << dataArray.size() << " ";
			// cout << tPB.x << " " << tPB.y << " " << tPB.z << endl;
			
			if (tPB.x == oldPB.x && tPB.y == oldPB.y && tPB.z == oldPB.z)
			{
				// Skip through
			}
			else
			{
				data_points.push_back(tPB);
				//inputAnalysis(tPB);
				
				//cout << data_points.size() << " ";
				//cout << tPB.x << " " << tPB.y << " " << tPB.z << endl;
				
				odx = dx;
				dx = (tPB.x - oldPB.x);
				ody = dy;
				dy = (tPB.y - oldPB.y);
				if ((odx == dx) && (ody == dy))
					step = dx;
				
				oldPB = tPB;
			}
   		}
	}
   	inData.close();
	
	inputAnalysis();
	
	if (compare_variables.translate) 
	{	
		if (compare_variables.centre_image) 
		{
			calcCOM(); // Calculate Centre of Mass
		}
		else 
		{
			translate(compare_variables.translate_x,compare_variables.translate_y);
		}
	}
	calculateDataArray();
	
	// Save Centred Image if need be
	if (compare_variables.translate && compare_variables.save_results)
	{
		if (compare_variables.centre_image)
		{
			saveImage("COM_" + name);
		}
		else 
		{
			saveImage("TRANS_" + name);
		}
	}
	
	return EXIT_SUCCESS;
}

void Information::inputAnalysis()
// After being passed a new dataArray this is used to work out totals, maximums, minimums, etc.
// No inputs or outputs as all data is held internally
{
	float st_save = step;
	int n_save = n;
	clearCounters();
	step = st_save;
	n = n_save;
	
	com_x.clear();
	com_y.clear();
	
	point tPB;
	tPB.x = tPB.y = tPB.z = 0;
	
	if (dataArray.size() > 0)
	{
		for (unsigned int i = 0; i < dataArray.size(); i++)
		{
			tPB.z = dataArray[i];
			tPB.x = (i%n)*step;
			tPB.y = (i/n)*step;
			/////////////////
			checkLimits(tPB);
		}
	}
	else if (data_points.size() > 0)
	{
		for (unsigned int i = 0; i < data_points.size(); i++)
		{
			tPB = data_points[i];
			// cout << i << " " << tPB.x << " " << tPB.y << " " << tPB.z << endl;
			/////////////////
			checkLimits(tPB);
		}
	}
	
	//calcM();
	// Calculates total width
	// Split up due to previous calculation errors
	m = (int) floor(((max_x-min_x+step)/step)+0.5);
	
	//calcN();
	// Calculates total height
	// Split up due to previous calculation errors
	n = (int) floor(((max_y-min_y+step)/step)+0.5); 
	
	// cout << m << " " << n << " " << m/step << " " << n/step << endl;
	// cout << max_x << " " << min_x << " " << max_y << " " << min_y << endl;
}

void Information::checkLimits(point tPB)
// Compare tPB to current limits
// Inputs(tPB) - point to be analysed
{
	// cout << tPB.x << " " << tPB.y << " " << tPB.z << endl;
	
	if (tPB.x > max_x)
	{
		max_x = tPB.x;
	}
	if (tPB.y > max_y)
	{
		max_y = tPB.y;
	}
	if (tPB.z > max_z)
	{
		max_z = tPB.z;
	}
	if (tPB.x < min_x)
	{
		min_x = tPB.x;
	}
	if (tPB.y < min_y)
	{
		min_y = tPB.y;
	}
	if (tPB.z < min_z)
	{
		min_z = tPB.z;
	}
	
	count++;
	x_total+=tPB.x;
	y_total+=tPB.y;
	z_total+=tPB.z;
	
	// This is for calculation of centre of mass of input //
	if (compare_variables.centre_image)
	{
		bool bCom_x = false;
		for (unsigned int i = 0; i < com_x.size(); i++)
		{
			if (com_x[i].x == tPB.x)
			{
				com_x[i].z+=tPB.z;
				bCom_x = true;
			}
		}
		
		if (!bCom_x)
		{
			com_x.push_back(tPB);
		}
		
		bool bCom_y = false;
		for (unsigned int i = 0; i < com_y.size(); i++)
		{
			if (com_y[i].y == tPB.y)
			{
				com_y[i].z+=tPB.z;
				bCom_y = true;
			}
		}
		
		if (!bCom_y)
		{
			com_y.push_back(tPB);
		}
	}
	///////////////////////////////////////////////////////
}

void Information::calcCOM()
// Calculates centre of mass, and shifts image if required
// Added 07/09/2010
{
	float x, y, mass, total;
	x = y = mass = total = 0;
	
	for (unsigned int i = 0; i < com_x.size(); i++)
	{
		mass += com_x[i].z;
		// Currently we have the value by pixels - commented out here gives the real value
		total += (((com_x[i].x-min_x)/step)*com_x[i].z);
		// total += (com_x[i].x*com_x[i].z);
	}
	
	// COM for X
	x = round(total/mass);
	mass = total = 0;
	
	for (unsigned int i = 0; i < com_y.size(); i++)
	{
		mass += com_y[i].z;
		// Currently we have the value by pixels - commented out here gives the real value
		total += (((com_y[i].y-min_y)/step)*com_y[i].z);
		// total += (com_y[i].y*com_y[i].z);
	}
	
	// COM for Y
	y = round(total/mass);
	
	// Calculate Centre of Grids
	float middle_x, middle_y;
	middle_x = round((max_x - min_x) / (2*step));
	middle_y = round((max_y - min_y) / (2*step));
	
	// Calculate Difference between COM and Centre of Grid
	
	float diff_x, diff_y;
	diff_x = middle_x - x;
	diff_y = middle_y - y;
	
	string xS;
	string yS;
	string middle_xS;
	string middle_yS;
	
	NumberToString(x,xS);
	NumberToString(y,yS);
	NumberToString(middle_x,middle_xS);
	NumberToString(middle_y,middle_yS);
	
#pragma omp critical
	{
		// Shift points around!
		// Points beyond original remit are moved to 0 values 
		output_content->push_back("");
		output_content->push_back("!!!!!!!!!!!!!!!!!! CENTRING IMAGE !!!!!!!!!!!!!!!!!!");
		output_content->push_back("Image COM X = " + xS + ", Y = " + yS + " Grid Centre X = " + middle_xS + ", Y = " + middle_yS);
		output_content->push_back("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	}
	
	translate(diff_x,diff_y);
}

void Information::translate(const float diff_x, const float diff_y)
// Translate data across array to new position
// Hopefully we'll be able to use this again
// Inputs - diff_x = shift in x axis
//		  - diff_y = shift in y axis
{
	string diff_xS;
	string diff_yS;
	
	NumberToString(diff_x,diff_xS);
	NumberToString(diff_y,diff_yS);
	
#pragma omp critical
	{
		output_content->push_back("!!!!!!!!!!!!!!!!!!!! NOTIFICATION !!!!!!!!!!!!!!!!!!!!");
		output_content->push_back("Shifting Image: X = +" + diff_xS + ", Y = +" + diff_yS);
		output_content->push_back("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	}
	
	//dataArray.clear(); // Have to remove all these points to start the array fresh
	
	if ((diff_x != 0) || (diff_y != 0))
	{	
		
		// cout << max_x << " " << min_x << " " << max_y << " " << min_y << endl;
		
		for (unsigned int i = 0; i < data_points.size(); i++)
		{
			data_points[i].x += diff_x;
			data_points[i].y += diff_y;
			
			// If beyond boundaries move to other side and set value to zero
			if (data_points[i].x > (max_x+(step/2)))
			{
				data_points[i].x -= (m*step);
				data_points[i].z = 0;
			}
			else if (data_points[i].x < (min_x-(step/2)))
			{
				data_points[i].x += (m*step);
				data_points[i].z = 0;
			}
			
			if (data_points[i].y > (max_y+(step/2)))
			{
				data_points[i].y -= (n*step);
				data_points[i].z = 0;
			}
			else if (data_points[i].y < (min_y-(step/2))) 
			{
				data_points[i].y += (n*step);
				data_points[i].z = 0;
			}
		}
	}
	
	// Analyse new data
	inputAnalysis();
	//calculateDataArray();
}

void Information::calculateDataArray()
// This will transfer the data from point format, i.e. just read, to an array, for quicker access
// Takes no inputs or outputs, just internal variables
{	
	//vector<point> backup_data_points = data_points;
	float position;
	point tempPB;
	// Initiate incase of misuse
	position = tempPB.x = tempPB.y = tempPB.z = 0;
	//
	unsigned int i, j, total_size;
	i = j = 0;
	total_size = data_points.size();
	
	while (i < total_size)
	{
		for (j = 0; j < data_points.size(); j++)
		{
			tempPB = data_points[j];
			position = ((((tempPB.y-min_y)/step)*n) + ((tempPB.x-min_x)/step));
			
			if (fabs(position-i) < ERRORBAR ) 
			{
				break;	
			}
		}
		dataArray.push_back(tempPB.z);
		data_points.erase(data_points.begin()+j);
		i++;
	}
	
	// Restore data points in case we need them in the future
	//data_points = backup_data_points;
}

void Information::scale(float scalar)
// Scale totals by amount requested. Best bet is max_z_first/max_z_second
// Inputs: float(scalar) - Value by which to multiply all intensities
{
	float newValue, oldValue;
	z_total = 0;
	max_z = 0;
	min_z = pow((float)rand(),2);
	
	for (int i = 0; i < count; i++)
	{
		oldValue = dataArray[i];
		newValue = oldValue*scalar;
		
		dataArray[i] = newValue;
		z_total+=newValue;
		
		if (newValue > max_z)
		{
			max_z = newValue;
		}
		
		if (newValue < min_z)
		{
			min_z = newValue;
		}
	}
}
