#ifndef CCLUSTER_H
#define CCLUSTER_H

/**
 02/06/2011
 - Commented out headers inherited from Utils.h
 **/

// #include <cstdlib>
// #include <string>
// #include <vector>
// #include <cmath>
#include <iomanip>
#include "Utils.h"
#include <fstream>
#include <iostream>

#define PI 3.141592653589793238462643383279502884197 //BIG PIE :)
#define DEG2RAD PI/180.0

// STRUCTURES //
struct atom
{
	std::string atom_type;
	double x;
	double y;
	double z;
	int coordination;
	int identity;
};

struct plane
{
	atom *patom_1;
	atom *patom_2;
	atom *patom_3;
	bool surface;
};

struct edge
{
	atom *patom_1;
	atom *patom_2;
};

struct vertex
{
	atom *patom_1;
};
///////////////

// CLUSTER CLASS DEFINITIONS //
class CCluster
{
public:
	CCluster();
	CCluster(std::string structure_filename);
	~CCluster() {;}
	
    //Cluster Physical properties
	std::vector<atom> get_atom_vec(){return atom_vec;}
	void set_atom_vec(std::vector<atom> vec){atom_vec = vec;}
	
	int get_nr_atoms(){return n_atoms;}
	int getAtomANum(){return n_atomsa;}
	int getAtomBNum(){return n_atomsb;}
	
	std::string getelementa(){return element_a;}
	std::string getelementb(){return element_b;}
	
	double get_energy(){return energy;}
	void set_energy(double val){energy = val;}
	
    //Transformations, Rotations and scaling (add scaling later) these are general transformation of points not cluster specific
	void rotate_x_axis(const double &radians);    
	void rotate_y_axis(const double &radians);    
	void rotate_z_axis(const double &radians);
	
	void translate(float x, float y, float z);   
	void shift_origin(atom *ptr);
	void shift_origin(atom &temp);
	void origin();
	void slice(const double &x, const double &y, const double &z);
	
	void scale(float x, float y, float z);
	void radial_Scale(float percentage);
	
    //cluster calculations
	void calc_surface_energy(); //*
	void calc_binding_energy(){binding_energy = energy/n_atoms;}
	double get_binding_energy(){calc_binding_energy(); return binding_energy;}
	
    //Some intersting methods
	
	double bond_length(const double &aX, const double &bX, const double &aY, const double &bY, const double &aZ, const double &bZ);
	std::vector<std::string> avg_bond_length(const std::string &a, const std::string &b, const float cutoff);
	double CalcMeanRadius();
	double CalcMeanRadius(const std::string &a);
	
    //More specfic movements
	void rotate_to_edge(int num); //This will rotate a selected edge to the xy plane *
	void rotate_to_plane(int num); //This will rotate the cluster on to an surface plane which one I dont know for Ico no big deal or planes equivalent *
	void place_vertex(int num); //this will place the cluster atom in the xy plane with the rest of the cluster above the plane *
	void print_surface_planes(); //*
	void print_edges(); //*
	void atom_type(); 
	
    //Sorting operators
	void place_furthest_atomic_distance_along_z_axis(); //Does exactly what it says on the tin *
	void place_two_atoms_in_line_z_axis(int atom1, int atom2); //*
	
	bool operator==(CCluster a); 
	bool same_atom(atom *aptr_1, atom *aptr_2); 
	friend std::ostream& operator<<(std::ostream& os, const CCluster& c);
	
private:
	double surface_area; //Does not take into account different metals
	
	std::string element_a;
	std::string element_b;
	
	int n_atoms;
	int n_atomsa;
	int n_atomsb;
	
	std::vector<plane> atom_plane;
	std::vector<plane> surface_planes;
	std::vector<edge>  m_edges;
	std::vector<atom> atom_vec;
	
	double energy;
	double binding_energy;
	
	float bond_scale_factor;
	
	void printatoms();
};
#endif
