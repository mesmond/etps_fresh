//***************************************************************************
//***************************************************************************
/*
 * File: mesmond-utils.h
 * Description: This file contains declarations of
 * 	several useful one-off functions. It includes the folowing:
 * 		min_2arg (double)
 * 		max_array (double)
 * 		isEqual (double)
*/
//***************************************************************************
//***************************************************************************
#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h> //exit(1)
#include <math.h>

using namespace std;

double min_2arg(double arg1, double arg2);
double max_2arg(double arg1, double arg2);
double max_array(double array[],int num_elements);

/*
 * Function: IsEqual()
 * Purpose: Compare floating point numbers for equality.
 *
 * Description: Adapted from
 * 	http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
 * 	This article explains that this particular approach will not behave
 * 	well with numbers very close to zero.
 */
bool isEqual(double A, double B, double maxRelativeError=1.0e-10);


#endif // INCLUDE_UTILS_H_
