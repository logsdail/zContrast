#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <cstdlib>
#include <string.h>
#include <vector>
#include <cmath>
#include <sstream>
#include "Structures.h"

#define tab '\t'
#define newline '\n'

template<class T>
void StringToNumber(const std::string &s, T &val)
{
	std::stringstream os;
	os << s;
	os >> val;
}

template<class T>
void NumberToString( T &val, std::string &s)
{
	std::stringstream os;
	os << val;
	os >> s;
}

struct DataPoint2D
{
	double x;
	double y;
};

struct vec3d
{
	double x;
	double y;
	double z;
	double mod;
};

void point_on_3D_line(const vec3d &b_vec, const vec3d &p2, vec3d &new_point, double nu);

void point_on_2D_line(const DataPoint2D &p1, const DataPoint2D &p2, DataPoint2D &new_point, double nu);

void find_white_space_end(const std::string &s, unsigned int &i);

void find_character_end(const std::string &s, unsigned int &i);

void Tokenize(const std::string& str,  std::vector<std::string>& tokens, const std::string& delimiters = " ");

bool cmpStr(const std::string &s1, const std::string &s2);

int randomNumber(int hi, int *idum);

float randomNumberF(float hi, int *idum);

float randomNumber(int *idum);

#endif
