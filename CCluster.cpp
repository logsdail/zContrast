#include "CCluster.h"

using namespace std;

ostream& operator<<(ostream& os, const CCluster& c)
// Inputs:  ostream(os) - Output stream (open)
//	    CCluster(c) - Cluster
// Outputs: ostream(os) - To file
{
	os << "  " << c.n_atoms << endl;
	os << setiosflags (ios::left|ios::fixed) << setprecision(8) << "   " << c.energy << endl;
	
	for(int i = 0; i < c.n_atoms; i++)
	{
		os << setiosflags ( ios::right|ios::fixed ) << setprecision(12) << c.atom_vec[i].atom_type << "  " 
		<< setw(16) << c.atom_vec[i].x * c.bond_scale_factor << "  "
		<< setw(16) << c.atom_vec[i].y * c.bond_scale_factor << "  " 
		<< setw(16) << c.atom_vec[i].z * c.bond_scale_factor;
		os << endl;
	}
	return os; 
}

bool CCluster::operator ==(CCluster a)
// Checks if cluster has rotated or scaled for tolerance
// Inputs: CCluster(a) - Cluster 
// Outputs: bool - Tolerance checking
{
	if(this->n_atoms != a.n_atoms) //Check to make sure cluster are of equal size first
	{
		return false;
	}
	
	this->origin();	// Move to Origin
	a.origin();	// Move to Origin
	
	for(int i = 0; i < this->n_atoms; i++) 	// Check for displacement will return false if rotation of one cluster has occured
	{
		double bond;
		bond = bond_length(this->atom_vec[i].x,a.atom_vec[i].x,this->atom_vec[i].y,a.atom_vec[i].y,this->atom_vec[i].z,a.atom_vec[i].z);
		
		if (bond >  0.5) return false;	// 0.5 is an arbitrary tolerance value. We are working in Angstrom
	}
	return true;
}

CCluster::CCluster()
// Initialise
{
	n_atoms = 0; // We can use this to tell if the object is empty
}

//CCluster::~CCluster()
//{
//	//delete &atom_vec;
//}

CCluster::CCluster(string structure_filename)
// Read Structure from file
// Input: string(structure_filename) - 
{
	string line;
	ifstream structure_file;
	atom temp;
	
	structure_file.open(structure_filename.c_str(), ios::in);
	if(!structure_file.is_open())
	{
		cout << "Structure file (" << structure_filename << ") does not exist" << endl;
		exit(1);
	}
	else
	{
		getline(structure_file, line); // Get number of atoms
		StringToNumber(line, n_atoms); 
		getline(structure_file, line); // Get energy
		if (line.size() != 0) StringToNumber(line, energy);

		while(!structure_file.eof())
		{
			vector<string> tokens;
			getline(structure_file, line);
			
			Tokenize(line, tokens, " \t");
			if (line.size() != 0)
			{
				temp.atom_type = tokens[0];
				StringToNumber(tokens[1], temp.x);
				StringToNumber(tokens[2], temp.y);
				StringToNumber(tokens[3], temp.z);
				atom_vec.push_back(temp);
			}
		}
		atom_type();
		bond_scale_factor = 1.0;
		structure_file.close();
	}
	//printatoms(); //Debugging
}

void CCluster::atom_type()
// Work out atom types
{
	string test;
	n_atomsb = n_atomsa = 0;        // Set number of atoms of each type
	element_a = test = atom_vec[0].atom_type;       // Set initial atom type
	
	for(int i = 0; i < n_atoms; i++)        // Work out atom types. Designed for bimetallic system, but will just do all atoms of first type and then all other atom types
		if(test == atom_vec[i].atom_type) n_atomsa++;
		else element_b = atom_vec[i].atom_type;
	
	n_atomsb = n_atoms - n_atomsa;
}

bool CCluster::same_atom(atom* aptr_1, atom* aptr_2)
// Check if two atoms are the same
// Input: atom(aptr_1) - Pointer to atom 1
//	  atom(aptr_2) - Pointer to atom 2
// Output: boolean
{
	if(aptr_1->atom_type != aptr_2->atom_type)
		return false; 
	
	if(aptr_1->identity != aptr_2->identity)
		return false;
	
	if(aptr_1->x != aptr_2->x)
		return false;
	
	if(aptr_1->y != aptr_2->y)
		return false;
	
	if(aptr_1->z != aptr_2->z)
		return false;
	
	return true;
}

double CCluster::bond_length(const double &aX, const double &bX, const double &aY, const double &bY, const double &aZ, const double &bZ)
// Calculate distance between two atoms
// Inputs: double(aX) - atom A, X coordinate
//	   double(bX) - atom B, X coordinate
//	   double(aY) - atom A, Y coordinate
//	   double(bY) - atom B, Y coordinate
//	   double(aZ) - atom A, Z coordinate
//	   double(bZ) - atom B, Z coordinate
// Outputs: double - atom displacements
{
	double bond, x_displacement, y_displacement, z_displacement;
	x_displacement = pow(aX-bX,2); // Square differences
	y_displacement = pow(aY-bY,2);
	z_displacement = pow(aZ-bZ,2);
	
	bond = sqrt(x_displacement+y_displacement+z_displacement);	// And square root total. Simple geometric algebra.
	
	return bond;
}

vector<string> CCluster::avg_bond_length(const string &a, const string &b, const float cutoff)
// Gives average bond length data
// Inputs: string(a) - atom A type
//	   string(b) - atom B type
//	   float(cutoff) - maximum distance we will pay attention to
// Outputs: vector<string> - List of strings with results to be displayed
{
	double bond , max_bond, min_bond, total;      
	min_bond = cutoff;
	bond = max_bond = total = 0.0;
	int count = 0;
	
	for(int first_atom = 0; first_atom < n_atoms; first_atom++)
		for(int second_atom = first_atom+1; second_atom < n_atoms; second_atom++)
			if( (atom_vec[first_atom].atom_type.compare(a)+atom_vec[second_atom].atom_type.compare(b) == 0) ||
			   (atom_vec[second_atom].atom_type.compare(b)+atom_vec[first_atom].atom_type.compare(a) == 0))
			{
				bond = bond_length(atom_vec[first_atom].x,atom_vec[second_atom].x,atom_vec[first_atom].y,atom_vec[second_atom].y,atom_vec[first_atom].z,atom_vec[second_atom].z);
				
				if (bond < cutoff)
				{
					total += bond;
					if (bond < min_bond) min_bond = bond;       
					if (bond > max_bond) max_bond = bond;
					count++;
				}
			} 
	
	// std::cout << "Atom A: " << a << std::endl;
	// std::cout << "Atom B: " << b << std::endl;
	// std::cout << "Bond Length Cutoff: " << cutoff << std::endl;
	// std::cout << "Number of A-B bonds: " << count << std::endl;
	// std::cout << "Average Bond Length: " << total/count << std::endl;
	// std::cout << "Minimum Bond Length: " << min_bond << std::endl;
	// std::cout << "Maximum Bond Length: " << max_bond << std::endl;
	
	string temp;
	vector<string> out;
	bond = total/count;
	
	temp = "Atom A: "; temp += a; out.push_back(temp);
	temp = "Atom B: "; temp += b; out.push_back(temp);
	NumberToString(cutoff,temp); temp = "Bond Length Cutoff: " + temp; out.push_back(temp);
	NumberToString(count,temp); temp = "Number of A-B Bonds: " + temp; out.push_back(temp);
	NumberToString(bond,temp); temp = "Average Bond Length: " + temp; out.push_back(temp);
	NumberToString(min_bond,temp); temp = "Minimum Bond Length: " + temp; out.push_back(temp);
	NumberToString(max_bond,temp); temp = "Maximum Bond Length: " + temp; out.push_back(temp);
	
	return out;
}

void CCluster::rotate_x_axis(const double &radians)
// Rotate cluster coordinates around x axis
// Input double(radians) - Degrees to rotate
{
	double temp_y, temp_z;
	for(int i = 0; i < n_atoms; i++)
	{
		temp_y = atom_vec[i].y * cos(radians) + atom_vec[i].z * sin(radians);
		temp_z = atom_vec[i].z * cos(radians) - atom_vec[i].y * sin(radians);
		atom_vec[i].y = temp_y;
		atom_vec[i].z = temp_z;
	}
}

void CCluster::rotate_y_axis(const double &radians)
// Rotate cluster coordinates around y axis
// Input double(radians) - Degrees to rotate
{
	double temp_x, temp_z;
	for(int i = 0; i < n_atoms; i++)
	{
		temp_x = atom_vec[i].x * cos(radians) - atom_vec[i].z * sin(radians);
		temp_z = atom_vec[i].x * sin(radians) + atom_vec[i].z * cos(radians);
		atom_vec[i].x = temp_x;
		atom_vec[i].z = temp_z;
	}
}

void CCluster::rotate_z_axis(const double &radians)
// Rotate cluster coordinates around z axis
// Input double(radians) - Degrees to rotate
{
	double temp_y, temp_x;
	for(int i = 0; i < n_atoms; i++)
	{
		temp_x = atom_vec[i].x * cos(radians) + atom_vec[i].y * sin(radians);
		temp_y = atom_vec[i].y * cos(radians) - atom_vec[i].x * sin(radians);
		atom_vec[i].y = temp_y;
		atom_vec[i].x = temp_x;
	}
}

void CCluster::shift_origin(atom *ptr)
// Shift origin to selected atom
// Input atom(ptr) - Pointer to atom to be used as origin
{
	origin();
	translate(ptr->x,ptr->y,ptr->z);
}

void CCluster::translate(float x, float y, float z)
// Translate atoms to new point 
// Input float(x) - X displacement
//	 float(y) - Y displacement
//	 float(z) - Z displacement
{
	for(int i = 0; i < n_atoms; i++)
	{
		atom_vec[i].x += x;
		atom_vec[i].y += y;
		atom_vec[i].z += z;
	}
}

void CCluster::shift_origin(atom &temp)
// Shift origin to selected atom
// Input atom(temp) - Atom to be referenced
{
	origin();
	translate(temp.x,temp.y,temp.z); // Get translate to do the work here
}

void CCluster::slice(const double &x, const double &y, const double &z)
// Slice molecule through give vector
// Input double(x) - x coordinate
//	 double(y) - y coordinate
//	 double(z) - z coordinate
{
	origin();
	
	std::vector<atom> SlicedModel;
	
	for(int i = 0; i < n_atoms; i++)
		if((atom_vec[i].x < x) &&  (atom_vec[i].y) < y && (atom_vec[i].z < z)){}
		else SlicedModel.push_back(atom_vec[i]);
	
	atom_vec = SlicedModel;
	n_atoms = SlicedModel.size();
	atom_type();
}

void CCluster::origin()
// Centralise cluster
{
	double avg_x, avg_y, avg_z;
	avg_x = avg_y = avg_z = 0;
	
	for( int i = 0; i < n_atoms; i++)
	{
		avg_x += atom_vec[i].x;
		avg_y += atom_vec[i].y;
		avg_z += atom_vec[i].z;
	}
	
	avg_x /= n_atoms;
	avg_y /= n_atoms;
	avg_z /= n_atoms;
	
	for( int i = 0; i < n_atoms; i++)
	{
		atom_vec[i].x -= avg_x;
		atom_vec[i].y -= avg_y;
		atom_vec[i].z -= avg_z;
	}
} 

double CCluster::CalcMeanRadius()
// Work out mean radius for whole cluster of certain type
// Output: double - average distance from origin to atoms
{
	double meanRadius;
	double radius_sum = 0.0;
	
	//Well before we do anything we best make sure the cluster COM is numerically found
	origin();
	
	for(int i = 0; i < n_atoms; i++)
		radius_sum += bond_length(atom_vec[i].x,0,atom_vec[i].y,0,atom_vec[i].z,0);
	
	meanRadius = (radius_sum/n_atoms);
	return meanRadius;
}

double CCluster::CalcMeanRadius(const string &a)
// Work out mean radius for atoms of certain type
// Input: string(a) - atom type
// Output: double - average distance from origin to atoms
{
	double meanRadius;
	double radius_sum = 0.0;
	int Na = 0;
   	
	//Well before we do anything we best make sure the cluster COM is numerically found
	origin();
	
	for(int i = 0; i < n_atoms; i++)
	{
		if( a == atom_vec[i].atom_type )
		{
			radius_sum += bond_length(atom_vec[i].x,0,atom_vec[i].y,0,atom_vec[i].z,0);
			Na++;      
		}
	}
	
	meanRadius = (radius_sum/Na);
	return meanRadius;   
}

void CCluster::scale(float x, float y, float z)
// Scale atom coordinates
// Input: float(x) - x scalar
//	  float(y) - y scalar
//	  float(z) - z scalar
{
	for(int currentAtom = 0; currentAtom < n_atoms; currentAtom++)
	{
		atom_vec[currentAtom].x *= x;
		atom_vec[currentAtom].y *= y;
		atom_vec[currentAtom].z *= z;
	}
}


void CCluster::radial_Scale(float percentage)
// Scale atom coordinates by percentage
// Input: float(percentage) - says on the tin
{
	double nu = 100 + percentage;
	vec3d origin,new_point,atom;
	origin.x = 0;
	origin.y = 0;
	origin.z = 0;
	
	for(int i = 0; i < n_atoms; i++)
	{
		atom.x = atom_vec[i].x;
		atom.y = atom_vec[i].y;
		atom.z = atom_vec[i].z;
		point_on_3D_line(origin,atom, new_point,nu);
		atom_vec[i].x = new_point.x;
		atom_vec[i].y = new_point.y;
		atom_vec[i].z = new_point.z;
	}
	
}

void CCluster::calc_surface_energy()
{
	double x,y,z,length,length1;
	for(int i = 0; i < n_atoms-2; i++)
	{
		for(int j = i+1; j < n_atoms-1; j++)
		{
			x = atom_vec[j].x - atom_vec[i].x;
			y = atom_vec[j].y - atom_vec[i].y;
			z = atom_vec[j].z - atom_vec[i].z;
			length = sqrt( x*x + y*y + z*z);
			
			if(length > 0 && length < 3.8) //REALLY NEED TO SET THIS TO A VARIABLE and not just 3.8
			{
				for(int k = j+1; k < n_atoms; k++)
				{
					x = atom_vec[k].x - atom_vec[i].x;
					y = atom_vec[k].y - atom_vec[i].y;
					z = atom_vec[k].z - atom_vec[i].z;
					length = sqrt( x*x + y*y + z*z);
					
					x = atom_vec[j].x - atom_vec[k].x;
					y = atom_vec[j].y - atom_vec[k].y;
					z = atom_vec[j].z - atom_vec[k].z;
					length1 = sqrt( x*x + y*y + z*z);
					
					//cout << "n_atoms = " << n_atoms << " i = " << i << " J = " << j << " K = " << k << endl;
					//cout << "Size of plane vector = " << atom_plane.size() << endl;
					
					if(length > 0 && length < 3.8 && length1 > 0 && length1 < 3.8)
					{
						//Then We have a candiate plane?
						plane temp_plane;
						temp_plane.patom_1 = &atom_vec[i];
						temp_plane.patom_2 = &atom_vec[j];
						temp_plane.patom_3 = &atom_vec[k];
						temp_plane.surface = false;
						atom_plane.push_back(temp_plane);
					}
					
				}
			}
		}
	}
	
	
	
	
	//===========================
	//All of the Above finds all nearest neighbour planes and adds them to a vector of planes
	//Now simple go through each plane and deduce if it is a surface or interior plane.
	//
	// CURRENT PROBLEMS WITH CODE BELOW THE ATOMVEC and PLANE data types dont comunicate
	//===========================
	
	double a,b, angleA, angleB, angleC, temp;
	angleA = angleB = angleC = 0.0;
	double x1, y1;
	atom cofm;
	
	for(unsigned int i = 0; i < atom_plane.size(); i++)
	{
		shift_origin(atom_plane[i].patom_1); //Shift the origin to atom on of plane 
		//Now work out the angle between the origin and atom2 and bring it in to the x axis using tangent
		a = atom_plane[i].patom_2->y * atom_plane[i].patom_2->y;
		b = atom_plane[i].patom_2->x * atom_plane[i].patom_2->x;
		a = sqrt(a);
		b = sqrt(b);
		
		temp = a/b;
		
		if( atom_plane[i].patom_2->x > 0.0 && atom_plane[i].patom_2->y > 0.0 )
		{
			//cout << "No Change (+x,+y)" << endl;
			angleA = atan(temp);
		}
		if( atom_plane[i].patom_2->x > 0.0 && atom_plane[i].patom_2->y < 0.0 )
		{
			//cout << "-1 * angle(+x,-y)" << endl;
			angleA = atan(temp);
			angleA = -1 * angleA; 
		}
		if( atom_plane[i].patom_2->x < 0.0 && atom_plane[i].patom_2->y < 0.0 )
		{
			//cout << "No Change (-x,-y)" << endl;
			angleA = atan(temp);	   
		}
		if( atom_plane[i].patom_2->x < 0.0 && atom_plane[i].patom_2->y > 0.0 )
		{
			//cout << "-1 * angle (-x,+y)" << endl;
			angleA = atan(temp);
			angleA = -1 * angleA; 
		}
		if( atom_plane[i].patom_2->y == 0.0 )
		{
			//cout << "Atom is 0 y z rotation" << endl;
			angleA = 0;
		}
		
		//cout << "----------------" << endl;
		//cout << "Before" << endl;
		//cout << atom_plane[i].patom_2->z << "\t" << atom_plane[i].patom_2->y << endl;
		
		rotate_z_axis(angleA); //Rotate system around
		
		//cout << "After" << endl;
		//cout << atom_plane[i].patom_2->z << "\t" << atom_plane[i].patom_2->y << endl;
		//cout << "----------------" << endl;
		//cout << endl;
		//Now remove the z- componet of atom 2 in plane i
		
		a = atom_plane[i].patom_2->z * atom_plane[i].patom_2->z;
		b = atom_plane[i].patom_2->x * atom_plane[i].patom_2->x;
		a = sqrt(a);
		b = sqrt(b);
		
		temp = a/b;
		
		if( atom_plane[i].patom_2->x > 0 && atom_plane[i].patom_2->z > 0 )
		{
			//cout << "360 - angleB (x+,z+)" << endl;
			angleB = atan(temp);
			angleB = -1*angleB;
		}
		if( atom_plane[i].patom_2->x > 0 && atom_plane[i].patom_2->z < 0 )
		{
			//cout << "No change (+x, -z)" << endl;
			angleB = atan(temp);
			//angleB = 360*DEGTORAD - 180*DEGTORAD + angleB; 
		}
		if( atom_plane[i].patom_2->x < 0 && atom_plane[i].patom_2->z < 0 )
		{
			//cout << "360 - angleB (-x, -z)" << endl;
			angleB = atan(temp);
			angleB = -1*angleB; 
		}
		if( atom_plane[i].patom_2->x < 0 && atom_plane[i].patom_2->z > 0 )
		{
			//cout << "No Change (-x, +z)" << endl;
			angleB = atan(temp);
			//angleB = 360*DEGTORAD - (180*DEGTORAD + angleB); 
		}
		if( atom_plane[i].patom_2->z == 0 )
		{
			//cout << "ZERO atom2 y rotation" << endl; 
			angleB = 0;
		}
		//cout << "Before" << endl;
		//cout << atom_plane[i].patom_2->y << "\t and z " << atom_plane[i].patom_2->z << endl;  
		
		rotate_y_axis(angleB); // Rotate system
		
		//cout << "after" << endl;
		//cout << atom_plane[i].patom_2->y << "\t and z " << atom_plane[i].patom_2->z << endl;
		//cout << endl;
		
		
		/*
		 cout << "Atom 3 " << atom_plane[i].patom_3->x << "\t" << atom_plane[i].patom_3->y << "\t" << atom_plane[i].patom_3->z << endl;
		 cout << "Atom 3 id" << atom_plane[i].patom_3->identity << endl;
		 cout << "Atom 3 by vec " << atom_vec[5].x << "\t" << atom_vec[5].y << "\t" << atom_vec[5].z << "\t" << atom_vec[5].identity << endl;
		 cout << "-------------" << endl;
		 */
		
		//Finally bring the 3 atom of plane i z componet to zero so now 3 atoms lie in the xy plane
		
		a = atom_plane[i].patom_3->z * atom_plane[i].patom_3->z;
		b = atom_plane[i].patom_3->y * atom_plane[i].patom_3->y;
		a = sqrt(a);
		b = sqrt(b);
		
		temp = a/b;
		
		if( atom_plane[i].patom_3->y > 0 && atom_plane[i].patom_3->z > 0 )
		{
			//cout << "no change (+y,+z)" << endl;
			angleC = atan(temp);
			//angleC = 360*DEGTORAD - (180*DEGTORAD + angleC);
		}
		if( atom_plane[i].patom_3->y > 0 && atom_plane[i].patom_3->z < 0 )
		{
			//cout << " -1*angleC (+y,-z)" << endl;
			angleC = atan(temp);
			angleC = -1*angleC; 
		}
		if( atom_plane[i].patom_3->y < 0 && atom_plane[i].patom_3->z < 0 )
		{
			//cout << "no change (-y,-z)" << endl;
			angleC = atan(temp);
			//angleC = 360*DEGTORAD + angleC; 
		}
		if( atom_plane[i].patom_3->y < 0 && atom_plane[i].patom_3->z > 0 )
		{
			//cout << "-1 * angleC (-y,+z)" << endl;
			angleC = atan(temp);
			angleC = -1* angleC; 
		}
		if( atom_plane[i].patom_3->z == 0 )
		{
			//cout << "ZERO atom_3 x rotation" << endl;
			angleC = 0;
		}
		
		//cout << "Before" << endl;
		//cout << "y = " << atom_plane[i].patom_3->y << " z = " << atom_plane[i].patom_3->z << endl; 
		
		rotate_x_axis(angleC); //Preform the rotation
		//DO YOU REMEBER WHY YOU DID ALL THE ABOVE BEN, I THINK NOT
		
		//cout << "After" << endl;
		//cout << "y = " << atom_plane[i].patom_3->y << " z = " << atom_plane[i].patom_3->z << endl;
		//cout << endl; 
		
		
		//Ok plane is now orientated by to x,y plane any atoms lying above and below z axis make this plane interior
		/*
		 cout << "Atom 1 " << atom_plane[i].patom_1->x << "\t" << atom_plane[i].patom_1->y << "\t" << atom_plane[i].patom_1->z << endl;
		 cout << "Atom 2 " << atom_plane[i].patom_2->x << "\t" << atom_plane[i].patom_2->y << "\t" << atom_plane[i].patom_2->z << endl;
		 cout << "Atom 3 " << atom_plane[i].patom_3->x << "\t" << atom_plane[i].patom_3->y << "\t" << atom_plane[i].patom_3->z << endl;
		 cout << "Atom 3 id" << atom_plane[i].patom_3->identity << endl;
		 cout << "Atom 3 by vec " << atom_vec[5].x << "\t" << atom_vec[5].y << "\t" << atom_vec[5].z << "\t" << atom_vec[5].identity << endl;
		 */
		//This make a new point which lies on the center of mass of the plane described by 3 atoms
		
		x1 = atom_plane[i].patom_1->x + atom_plane[i].patom_2->x + atom_plane[i].patom_3->x;
		y1 = atom_plane[i].patom_1->y + atom_plane[i].patom_2->y + atom_plane[i].patom_3->y;
		
		cofm.x = x1/3.0;
		cofm.y = y1/3.0;
		cofm.z = 0.0;
		
		shift_origin(cofm);
		
		float near_zero = 1e-12; //anything less than this is smaller than the orignal coordination file input.
		
		bool bIsNegative = false;
		bool bIsPositive = false;
		
		for(int j = 0; j < n_atoms; j++)
		{
			float distancexy; //Dilerberate lose of precision its a cutoff point if it 3.8 dont need all them zeros
			float distancexyz, radius_cofm;
			radius_cofm = 1.615;//THIS IS NOT IDEAL REALLY SHOULD BE CALCULATED FROM Atom Componenets SEE LAB BOOK
			
			distancexy = atom_vec[j].x * atom_vec[j].x + atom_vec[j].y * atom_vec[j].y;
			distancexy = sqrt(distancexy);
			distancexyz = atom_vec[j].x * atom_vec[j].x +  atom_vec[j].y * atom_vec[j].y + atom_vec[j].z * atom_vec[j].z;
			
			if( distancexy < radius_cofm )//then it lies inside of plane radius cutoff
			{
				if( atom_vec[j].z < -near_zero || atom_vec[j].z > near_zero ) //Then there is an atom lying on onside of the plane
				{
					if( atom_vec[j].z < -near_zero ) //is it the negative side?
					{
						//cout << "Plane " << i << " has negative side" << endl;
						bIsNegative = true;
					}
					if( atom_vec[j].z > near_zero ) // is it postive side ?
					{
						//cout << "Plane " << i << " has postive side" << endl;
						bIsPositive = true;
					}
					if( bIsPositive && bIsNegative ) //Interior plane
					{
						atom_plane[i].surface = false;
						break;
					}
					else
					{
						atom_plane[i].surface = true;
					}
				}
			}
		}
	}
	int surfaces = 0;
	double surfacearea = 0;
	for(unsigned int i = 0; i < atom_plane.size(); i++)
	{
		if(atom_plane[i].surface == true) //If the plane was found to be a surface plane
		{
			surface_planes.push_back(atom_plane[i]);
			double a1, b1, c1, C;
			a1=b1=c1=C=0.0;
			//Area a*bsin(C)/2  area of a general triangle
			a1 = (atom_plane[i].patom_1->x - atom_plane[i].patom_2->x)*(atom_plane[i].patom_1->x - atom_plane[i].patom_2->x) + 
			(atom_plane[i].patom_1->y - atom_plane[i].patom_2->y)*(atom_plane[i].patom_1->y - atom_plane[i].patom_2->y) +
			(atom_plane[i].patom_1->z - atom_plane[i].patom_2->z)*(atom_plane[i].patom_1->z - atom_plane[i].patom_2->z);
			a1 = sqrt(a1);
			b1 = (atom_plane[i].patom_1->x - atom_plane[i].patom_3->x)*(atom_plane[i].patom_1->x - atom_plane[i].patom_3->x) + 
			(atom_plane[i].patom_1->y - atom_plane[i].patom_3->y)*(atom_plane[i].patom_1->y - atom_plane[i].patom_3->y) +
			(atom_plane[i].patom_1->z - atom_plane[i].patom_3->z)*(atom_plane[i].patom_1->z - atom_plane[i].patom_3->z);
			b1 = sqrt(b1);
			c1 = (atom_plane[i].patom_3->x - atom_plane[i].patom_2->x)*(atom_plane[i].patom_3->x - atom_plane[i].patom_2->x) + 
			(atom_plane[i].patom_3->y - atom_plane[i].patom_2->y)*(atom_plane[i].patom_3->y - atom_plane[i].patom_2->y) +
			(atom_plane[i].patom_3->z - atom_plane[i].patom_2->z)*(atom_plane[i].patom_3->z - atom_plane[i].patom_2->z);
			c1 = sqrt(c1);
			
			C = acos((a1*a1 + b1*b1 - c1*c1)/(2*a1*b1)); //law of cosines finds angle C which 
			
			//cout << "a1 = " << a1 << endl;
			//cout << "b1 = " << b1 << endl;
			//cout << "c1 = " << c1 << endl;
			//cout << "C = " << C << endl;
			surfacearea += a1 * b1 * sin(C)/ 2.0; //Area of that triangle;
			surfaces++;
		}
	}
	//cout << "Number of surface planes = " << surfaces << " out of " << atom_plane.size() << " planes " << structure << endl;
	//cout << "Surface area = " << surfacearea << " Angstroms Squared " << endl;
	surface_area = surfacearea;
}

void CCluster::print_surface_planes()
{
	cout << "Plane" << "\t" << "atom 1" << "\t" << "atom 2" << "\t" << "atom 3" << endl;
	for(unsigned int i = 0; i < surface_planes.size(); i++)
	{
		cout << i << "\t" << surface_planes[i].patom_1->identity << "\t" << surface_planes[i].patom_2->identity << "\t" 
		<< surface_planes[i].patom_3->identity
		<< endl;
	}  
}

void CCluster::place_two_atoms_in_line_z_axis(int atom1, int atom2)
{
	atom temp_atom = atom_vec[atom1];
	shift_origin(temp_atom); //move origin to atom1
	
	double a,b,angleA, angleB, temp;
	
	a = atom_vec[atom2].y * atom_vec[atom2].y;
	b = atom_vec[atom2].x * atom_vec[atom2].x;
	
	a = sqrt(a);
	b = sqrt(b);
	
	temp = a/b;
	
	if(atom_vec[atom2].x > 0.0 && atom_vec[atom2].y > 0.0) //We are in +x +y space
	{
		angleA = atan(temp); 
	}
	if(atom_vec[atom2].x > 0.0 && atom_vec[atom2].y < 0.0) //We are in +x -y space
	{
		angleA = atan(temp);
		angleA = -1.0 * angleA; 
	}
	if(atom_vec[atom2].x < 0.0 && atom_vec[atom2].y < 0.0) //We are in -x -y space
	{
		angleA = atan(temp); 
	}
	if(atom_vec[atom2].x < 0.0 && atom_vec[atom2].y > 0.0) //We are in -x +y space
	{
		angleA = atan(temp); 
		angleA = -1.0 * angleA;
	}
	if(atom_vec[atom2].y == 0.0)
	{
		angleA = 0;
	}
	
	rotate_z_axis(angleA);
	
	a = atom_vec[atom2].z * atom_vec[atom2].z;
	b = atom_vec[atom2].x * atom_vec[atom2].x;
	
	a = sqrt(a);
	b = sqrt(b);
	
	temp = a/b;
	
	if(atom_vec[atom2].x > 0.0 && atom_vec[atom2].z > 0.0) //We are in +x +y space
	{
		angleB = atan(temp); 
	}
	if(atom_vec[atom2].x > 0.0 && atom_vec[atom2].z < 0.0) //We are in +x -y space
	{
		angleB = atan(temp);
		angleB = -1.0 * angleB; 
	}
	if(atom_vec[atom2].x < 0.0 && atom_vec[atom2].z < 0.0) //We are in -x -y space
	{
		angleB = atan(temp); 
	}
	if(atom_vec[atom2].x < 0.0 && atom_vec[atom2].z > 0.0) //We are in -x +y space
	{
		angleB = atan(temp); 
		angleB = -1.0 * angleB;
	}
	if(atom_vec[atom2].z == 0.0)
	{
		angleB = 0;
	}
	
	rotate_y_axis(angleB); // Rotate system
	
	origin();
	
	rotate_y_axis(90*DEG2RAD);
}


void CCluster::place_furthest_atomic_distance_along_z_axis()
{
	int atom1, atom2;
	atom1 = atom2 = 0;
	double distance, temp_distance;
	distance = 0.0;
	
	for( int i = 0; i < n_atoms; i++)
		for( int j = i+1; j < n_atoms; j++)
		{
			temp_distance = bond_length(atom_vec[i].x,atom_vec[j].x,atom_vec[i].y,atom_vec[j].y,atom_vec[i].z,atom_vec[j].z);
			if(temp_distance > distance)
			{
				distance = temp_distance;
				atom1 = i;
				atom2 = j;
			}
		}    
	//cout << atom1 << "\t" << atom2 << "\t" << distance << endl;
	
	place_two_atoms_in_line_z_axis(atom1,atom2);
}

void CCluster::rotate_to_plane(int num)
{  
	plane p;
	double a ,b , angleA, angleB, angleC, temp;
	angleA = angleB = angleC = 0.0;
	double x1, y1;
	atom cofm;
	p = surface_planes[num];
	//Rotate to a surface plane
	shift_origin(p.patom_1); //Shift the origin to atom on of plane i
	
	//Now work out the angle between the origin and atom2 and bring it in to the x axis using tangent
	a = p.patom_2->y * p.patom_2->y;
	b = p.patom_2->x * p.patom_2->x;
	a = sqrt(a);
	b = sqrt(b);
	
	temp = a/b;
	
	if( p.patom_2->x > 0.0 && p.patom_2->y > 0.0 )
	{
		//cout << "No Change (+x,+y)" << endl;
		angleA = atan(temp);
	}
	if( p.patom_2->x > 0.0 && p.patom_2->y < 0.0 )
	{
		//cout << "-1 * angle(+x,-y)" << endl;
		angleA = atan(temp);
		angleA = -1 * angleA; 
	}
	if( p.patom_2->x < 0.0 && p.patom_2->y < 0.0 )
	{
		//cout << "No Change (-x,-y)" << endl;
		angleA = atan(temp);	   
	}
	if( p.patom_2->x < 0.0 && p.patom_2->y > 0.0 )
	{
		//cout << "-1 * angle (-x,+y)" << endl;
		angleA = atan(temp);
		angleA = -1 * angleA; 
	}
	if( p.patom_2->y == 0.0 )
	{
		//cout << "Atom is 0 y z rotation" << endl;
		angleA = 0;
	}
	
	//cout << "----------------" << endl;
	//cout << "Before" << endl;
	//cout << p.patom_2->z << "\t" << p.patom_2->y << endl;
	
	rotate_z_axis(angleA); //Rotate system around
	
	//cout << "After" << endl;
	//cout << p.patom_2->z << "\t" << p.patom_2->y << endl;
	//cout << "----------------" << endl;
	//cout << endl;
	//Now remove the z- componet of atom 2 in plane i
	
	a = p.patom_2->z * p.patom_2->z;
	b = p.patom_2->x * p.patom_2->x;
	a = sqrt(a);
	b = sqrt(b);
	
	temp = a/b;
	
	if( p.patom_2->x > 0 && p.patom_2->z > 0 )
	{
		//cout << "360 - angleB (x+,z+)" << endl;
		angleB = atan(temp);
		angleB = -1*angleB;
	}
	if( p.patom_2->x > 0 && p.patom_2->z < 0 )
	{
		//cout << "No change (+x, -z)" << endl;
		angleB = atan(temp);
		//angleB = 360*DEGTORAD - 180*DEGTORAD + angleB; 
	}
	if( p.patom_2->x < 0 && p.patom_2->z < 0 )
	{
		//cout << "360 - angleB (-x, -z)" << endl;
		angleB = atan(temp);
		angleB = -1*angleB; 
	}
	if( p.patom_2->x < 0 && p.patom_2->z > 0 )
	{
		//cout << "No Change (-x, +z)" << endl;
		angleB = atan(temp);
		//angleB = 360*DEGTORAD - (180*DEGTORAD + angleB); 
	}
	if( p.patom_2->z == 0 )
	{
		//cout << "ZERO atom2 y rotation" << endl; 
		angleB = 0;
	}
	//cout << "Before" << endl;
	//cout << p.patom_2->y << "\t and z " << p.patom_2->z << endl;  
	
	rotate_y_axis(angleB); // Rotate system
	
	//cout << "after" << endl;
	//cout << p.patom_2->y << "\t and z " << p.patom_2->z << endl;
	//cout << endl;
	
	
	/*
	 cout << "Atom 3 " << p.patom_3->x << "\t" << p.patom_3->y << "\t" << p.patom_3->z << endl;
	 cout << "Atom 3 id" << p.patom_3->identity << endl;
	 cout << "Atom 3 by vec " << atom_vec[5].x << "\t" << atom_vec[5].y << "\t" << atom_vec[5].z << "\t" << atom_vec[5].identity << endl;
	 cout << "-------------" << endl;
	 */
	
	//Finally bring the 3 atom of plane i z componet to zero so now 3 atoms lie in the xy plane
	
	a = p.patom_3->z * p.patom_3->z;
	b = p.patom_3->y * p.patom_3->y;
	a = sqrt(a);
	b = sqrt(b);
	
	temp = a/b;
	
	if( p.patom_3->y > 0 && p.patom_3->z > 0 )
	{
		//cout << "no change (+y,+z)" << endl;
		angleC = atan(temp);
		//angleC = 360*DEGTORAD - (180*DEGTORAD + angleC);
	}
	if( p.patom_3->y > 0 && p.patom_3->z < 0 )
	{
		//cout << " -1*angleC (+y,-z)" << endl;
		angleC = atan(temp);
		angleC = -1*angleC; 
	}
	if( p.patom_3->y < 0 && p.patom_3->z < 0 )
	{
		//cout << "no change (-y,-z)" << endl;
		angleC = atan(temp);
		//angleC = 360*DEGTORAD + angleC; 
	}
	if( p.patom_3->y < 0 && p.patom_3->z > 0 )
	{
		//cout << "-1 * angleC (-y,+z)" << endl;
		angleC = atan(temp);
		angleC = -1* angleC; 
	}
	if( p.patom_3->z == 0 )
	{
		//cout << "ZERO atom_3 x rotation" << endl;
		angleC = 0;
	}
	
	//cout << "Before" << endl;
	//cout << "y = " << p.patom_3->y << " z = " << p.patom_3->z << endl; 
	
	rotate_x_axis(angleC); //Preform the rotation
	//DO YOU REMEBER WHY YOU DID ALL THE ABOVE BEN, I THINK NOT
	
	//cout << "After" << endl;
	//cout << "y = " << p.patom_3->y << " z = " << p.patom_3->z << endl;
	//cout << endl; 
	
	
	//Ok plane is now orientated by to x,y plane any atoms lying above and below z axis make this plane interior
	/*
	 cout << "Atom 1 " << p.patom_1->x << "\t" << p.patom_1->y << "\t" << p.patom_1->z << endl;
	 cout << "Atom 2 " << p.patom_2->x << "\t" << p.patom_2->y << "\t" << p.patom_2->z << endl;
	 cout << "Atom 3 " << p.patom_3->x << "\t" << p.patom_3->y << "\t" << p.patom_3->z << endl;
	 cout << "Atom 3 id" << p.patom_3->identity << endl;
	 cout << "Atom 3 by vec " << atom_vec[5].x << "\t" << atom_vec[5].y << "\t" << atom_vec[5].z << "\t" << atom_vec[5].identity << endl;
	 */
	//This make a new point which lies on the center of mass of the plane described by 3 atoms
	
	x1 = p.patom_1->x + p.patom_2->x + p.patom_3->x;
	y1 = p.patom_1->y + p.patom_2->y + p.patom_3->y;
	
	//cofm.x = x1/3.0;
	//cofm.y = y1/3.0;
	//cofm.z = 0.0;
	
	//shift_origin(cofm);
	origin();
	rotate_y_axis(180*DEG2RAD);  //Not sure if this is sensible
	
}
/*  Not so bad as all surface planes can be identified edges between polygons can be identified, the actual edge (ie between to polyhedra vertecices)
 will be a liner combination of polygon edges in the same direction but for this purpose visually identifying an edge is all that is necessary.
 
 */
void CCluster::print_edges()
{
	vector<edge> ext_edges;
    //ext_edges = vector<edge>(surface_planes.size()*3); //3 atoms per face and all atoms are involved in edges but over estimation on number edges
	edge temp_edge;
	
	for(unsigned int i = 0; i < surface_planes.size(); i++)
	{
		
		temp_edge.patom_1 = surface_planes[i].patom_1;
		temp_edge.patom_2 = surface_planes[i].patom_2;
		ext_edges.push_back(temp_edge);
		
		temp_edge.patom_1 = surface_planes[i].patom_1;
		temp_edge.patom_2 = surface_planes[i].patom_3;
		ext_edges.push_back(temp_edge);
		
		temp_edge.patom_1 = surface_planes[i].patom_3;
		temp_edge.patom_2 = surface_planes[i].patom_2;
		ext_edges.push_back(temp_edge);
		
	}  
	
	cout << ext_edges.size() << endl;
	m_edges.push_back(ext_edges[0]);
	
	for(unsigned int i = 1; i < ext_edges.size(); i++)
	{
		temp_edge = ext_edges[i];
		for(unsigned int j = 0; j < m_edges.size(); j++)
		{
			if( same_atom(m_edges[j].patom_1, temp_edge.patom_1) && same_atom(m_edges[j].patom_2, temp_edge.patom_2))
			{
				
				cout << j << "TEST" << endl;
			}
			else
			{
				m_edges.push_back(temp_edge);      
			}
		}
		cout << "------------" << endl;
	} 
	cout << m_edges.size() << endl;
	for(unsigned int i = 0; i < m_edges.size(); i++)
	{
		cout << "Edge\t Atom_1\t Atom_2" << endl; 
		cout << i << "\t" << m_edges[i].patom_1->identity << "\t" << m_edges[i].patom_2->identity << endl; 
	}    
}

void CCluster::rotate_to_edge(int num)
{
	
}

/* This one is gonna be hard to generalise easy for ico but you have to identify
 Vertex which are on a surface plane and have a large number of of edges. Of course various polyhedra
 dont have to have same number of edges and also as an atom represent a vertex coordinate other atoms which
 lie in planes or on edges between to vertices of the polyhedra have to be excluded.
 */
void CCluster::place_vertex(int num)
{
	
}

void CCluster::printatoms()
{
	atom temp;
	for (int i = 0; i < n_atoms; i++)
	{
		temp = atom_vec[i];
		cout << temp.atom_type << " " << temp.x << " " << temp.y << " " << temp.z << endl;
	}
}
