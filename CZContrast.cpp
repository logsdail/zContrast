#include "CZContrast.h"

using namespace std;

void CZcontrast::init(int *s, vector<string> *o)
// Constructor for this intensity calculator
{
	// Initiate Random Numbers
	idum = s;
	output_content = o;
	//////////////////////////
	
	m_iNumY = m_iNumX = 0;
	m_fProbeSize = 0.0;
	m_fAlpha = 0.0;
	m_fExponent = m_fNoise = 0.0;
	
	// These are for soft coding of parameter componenets
	m_fAtomicNumber.clear();
	m_fElementName.clear();
	m_fAtomicRadii.clear();
}

CZcontrast::~CZcontrast()
// Deconstructor - deletes image
{
	delete m_pImage;
}

void CZcontrast::Image()
// Calculates image intensities
{
	std::vector<atom> temp;
	if(m_pCluster->get_nr_atoms() != 0)
	{
		m_pImage = new CContrastImage(m_iNumX, m_iNumY, m_fGridSize, output_content); //create a new image
		atom probe; //create an atom which we will use as a probe 
		
		probe.x = 0.5 * (m_iNumX * m_fGridSize); //this should give us the center point of the grid
		probe.y = 0.5 * (m_iNumY * m_fGridSize);
		probe.z = 0.0; //no 3rd Dimension
		
		m_pCluster->shift_origin(probe); //Now that the probe is central make its position the clusters center
		temp = m_pCluster->get_atom_vec(); //ok now copy positions of cluster into this atom vector	
		
		//new time to move the probe in left -> right direction and down after reseting it.
		for(int i = 0; i < m_iNumY; i++)
			for(int j = 0; j < m_iNumX; j++)
			{
				vector<double> atom_count;
				
				if (m_fElementName.size() == 0) 
				{
					atom_count.resize(1);
				}
				else 
				{
					atom_count.resize(m_fElementName.size());
				}
				
				double intensity = 0.0;
				
				//ok move probe to correct position 
				probe.x = j * m_fGridSize;
				probe.y = i * m_fGridSize;
				
				SearchForAtoms(temp,probe, atom_count); // Search for atoms
				intensity = DetectIntensity(atom_count); // Work out intensity at this grid point
				if (m_fNoise > 0)
				{
					intensity = addNoise(intensity);
				}
				// m_pImage->setIntensity(j,i,intensity); // Add to data
				m_pImage->setIntensity(i,j,intensity);	  // X- and Y- incorrectly defined above. Caused 90 degree rotation
			}
	}
}

void CZcontrast::setImageArea(float x, float y)
// Set overall image area by multiplying grid size by int
// Inputs: float(x) - Number of points on X axis
//	   float(y) - Number of points on y axis
{
	m_fPointsPerLength = 1.0/m_fGridSize; 
	m_iNumX = (int) ceil((x * m_fPointsPerLength) - 0.5);
	m_iNumY = (int) ceil((y * m_fPointsPerLength) - 0.5);
}

void CZcontrast::SearchForAtoms(vector<atom> vec, atom &probe, vector<double> &atom_count)
// Go through atoms and calculate their intensity
// Inputs: vector(atoms) - All atoms in cluster
//	   atom(probe) - Probe location
//	   vector(atom_count) - Saves atom intensity for each location to be multiplied by atom type
{
	// Go through atoms and look at intensity
	for(vector<atom>::iterator it = vec.begin(); it != vec.end(); it++)
	{
		// Work out distance between atom and probe
		double distance = ((*it).x - probe.x)*((*it).x - probe.x);
		distance += ((*it).y - probe.y)*((*it).y - probe.y);
		distance = sqrt(distance);
		
		// Next we check proximity to probe.
		// As we are dealing with small clusters I am going to remove this, but it could be softcoded in the future
		// as perhaps a cutoff distance. Currently 300Angrstrom is well within our image sizes
		
		bool flag = false;
		unsigned int i = 0;
		
		if (m_fElementName.size() == 0) 
		{
			atom_count[0] += Gaussian(distance,9999);
		}
		// Work out atom type and thus contribution to intensity
		else
		{
			while (!flag) 
			{
				if((*it).atom_type == m_fElementName[i])
				{
					atom_count[i] += Gaussian(distance,i);
					flag = !flag;
				}
				else
				{
					i++;
					if(i >= m_fElementName.size())
					{
#pragma omp critical
						{
							output_content->push_back("This type of element is not defined in the parameters:" + (*it).atom_type);						
							output_content->push_back("Elementsfile name must contain in the correct order:");
							output_content->push_back("(int) Number of Elements defined in parameter file");
							output_content->push_back("(str) Element");
							output_content->push_back("(int) Atomic Number");
							output_content->push_back("(float) Atomic Radii");
							output_content->push_back("(Last three repeated for however many elements documented)");
						}
						flag = !flag;
					}
				}
			}
		}
	}
}

double CZcontrast::Gaussian(const double &distance,const int &i)
// Calculate general intensity from Gaussian function
// Inputs double(distance) - Distance from point
// 	  int(i) - Atom number in cluster
// Outputs: double - Overall intensity without scaling
{
	double intensity = 0.0;
	double x = 0.0;
	double radiiTemp = 1.0;
	double radiiFactor = 0.0;
	
	x = distance*distance; //Square distance
	
	// If i = 9999 we have a are not using atomic radius
	if (i != 9999) 
	{
		radiiTemp = 1/m_fAtomicRadii[i];
	}
	
	radiiFactor = radiiTemp*radiiTemp; // Square radius over one
	
	intensity = exp(-m_fAlpha*x*radiiFactor); // e^alpha*(d/r)^2
	return intensity;
}

double CZcontrast::DetectIntensity(vector<double> atom_count)
// Calculate scaling for intensity
// Input vector(atom_count) - Each atoms gaussian contribution
// Outputs: Their collective intensity
{
	
	double intensity = 0;
	
	if (m_fAtomicNumber.size() == 0) 
	{
		intensity = atom_count[0];
	}
	else
	{
		// Sum intensity of each atom combined together at current point
		for (unsigned int i = 0; i < atom_count.size(); i++)
		{
			intensity += pow(m_fAtomicNumber[i], m_fExponent) * atom_count[i];
		}
	}
	
	return intensity;
}

double CZcontrast::addNoise(const double &a)
// Add noise to value, via random number generator
// Inputs: float a - starting value
// Outputs: float with random noise added
// Added 13/10/2010
{
	// double i = (m_fNoise - randomNumberF(m_fNoise*2));
	double i = randomNumberF(m_fNoise,idum);
	
	i += a;
	
	return i;
}

/**
 void CZcontrast::writeCrossSection(string filename)
 // Write cross section
 // Inputs: String(filename) -
 {
 ofstream output;
 output.open(filename.c_str(),ios::out);
 
 if(!output.is_open())
 cout << "Error Writing Contrast File Skipping CrossSection" << std::endl;
 else
 {
 // Find highest intensity to normalise data
 // float bestIntensity = 0.0;
 
 // for(size_t currentData = 0; currentData < m_vData.size(); currentData++)
 // {
 // 	if(m_vData[currentData].y > bestIntensity)
 // 	{
 // 		bestIntensity = m_vData[currentData].y;
 //	}
 // }
 
 //output in NanoMeters and Normalised by highest intensity.
 for(size_t currentData = 0; currentData < m_vData.size(); currentData++)
 {
 output << m_vData[currentData].x << "\t" << m_vData[currentData].y << endl;
 // output << m_vData[currentData].x << "\t" << m_vData[currentData].y / bestIntensity << endl;
 }
 }
 output.close();
 }
 **/

/**
 void CZcontrast::CrossSection(float x1, float y1, float x2, float y2)
 // Cross section calculator. Use can be softcoded?
 // Inputs float(x1) - Start point X coordinate
 //        float(y1) - Start point Y coordinate
 //        float(x2) - Finish Point X coordinate
 //        float(y2) - Finish point Y coordinate
 // Outputs: To m_vData for printing
 {
 m_vData.clear(); //just making sure no false crosssection data
 vector<atom> temp;
 data currentRecord;
 
 // PROBLEM IS HERE IN THIS METHOD WHEN USING IN CURRENT LOCATION
 // SEEMS TO HAVE STOPPED WORKING NOW ROTATION HAS BEEN INTRODUCED
 // NOT TO WORRY AS WE AREN"T USING THIS ANYWAY
 temp = m_pCluster->get_atom_vec();
 
 DataPoint2D p1, p2, new_point;
 double distance, nu;
 atom probe;
 
 probe.x = p1.x = x1; // starting point for probe
 probe.y = p1.y = y1; // starting point for probe
 probe.z = 0.0;
 p2.x = x2;
 p2.y = y2;
 
 distance  = (p2.x - p1.x)*(p2.x - p1.x);
 distance += (p2.y - p1.y)*(p2.y - p1.y);
 distance = sqrt(distance); // Distance between points
 int num_points = (int) (floor(distance)/m_fGridSize); // divide distance into segments
 nu = 1.0/num_points;
 
 for(int i = 0; i < num_points; i++)
 {
 vector<double> atom_count;
 
 if (m_fElementName.size() == 0) 
 atom_count.resize(1);
 else 
 atom_count.resize(m_fElementName.size());
 
 point_on_2D_line(p1,p2,new_point, (nu*i));
 probe.x = new_point.x; // Move probe along line
 probe.y = new_point.y;
 
 SearchForAtoms(temp, probe, atom_count); // Count atoms
 
 currentRecord.x = i*m_fGridSize;
 currentRecord.y = DetectIntensity(atom_count);
 
 m_vData.push_back(currentRecord); //Store information
 }
 }
 **/
