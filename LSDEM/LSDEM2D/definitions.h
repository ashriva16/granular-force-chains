/*
 * definitions.h
 *
 *  Created on: July 14, 2014
 *      Author: Reid
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

// eigen doesn't support vectorization for window$
// #ifdef _WIN32
// 	#define EIGEN_DONT_ALIGN_STATICALLY
// #endif

// C libraries
#include <float.h>		// DBL_MAX
#include <math.h> 		// sin, cos
#include <stdio.h>      // printf, fgets
#include <stdlib.h>     // atoi
#include <omp.h>		// OpenMP
#include <time.h>		// time
#include <mpi.h>		// MPI

// C++ libraries
#include <string>		// strings
using std::string;

#include <iostream>		// terminal output
#include <cmath>
using std::cout;
using std::endl;

#include <fstream>		// file io
using std::ifstream;
using std::getline;
using std::istringstream;

#include <algorithm>	// math stuff
using std::min;
using std::max;
using std::find;

#include <vector>  // standard vector
//#include <array>
using std::vector;
//using std::array;
 
#include <ctime>
using std::clock_t;

#include<tr1/array>
using std::tr1::array;

#include <iomanip>

#include <Eigen>
using Eigen::Vector2d;
using Eigen::Vector2i;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix2d;
using Eigen::Matrix;
typedef Matrix<double, 8, 1> Vector8d;
typedef Matrix<double, 6, 1> Vector6d;

 

#include <cstdlib>
using std::rand;


#endif /* DEFINITIONS_H_ */
