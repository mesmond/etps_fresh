//***************************************************************************
//***************************************************************************
/*
 * File: mesmond-utils.h
 * Description: This file contains declarations and definitions of
 * 	several useful one-off functions. It includes the folowing:
*/
//***************************************************************************
//***************************************************************************
#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h> //exit(1)

using namespace std;

double min_2arg(double arg1, double arg2)
{
	double min=arg1;
	if (arg2<min) min=arg2;
	return min;
}

double max_array(double array[],int num_elements)
{
	int i;
	double max=array[0];
	for (i=0;i<num_elements;++i)
	{
		if (array[i]>max)
		{
			max=array[i];
		}
	}
	return max;
}


/*
 * Function: IsEqual()
 * Purpose: Compare floating point numbers for equality.
 *
 * Description: Adapted from
 * 	http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
 * 	This article explains that this particular approach will not behave
 * 	well with numbers very close to zero.
 */
bool isEqual(double A, double B, double maxRelativeError=1.0e-10)
{
    if (A == B)
        return true;
    double relativeError;
    if (fabs(B) > fabs(A))
        relativeError = fabs((A - B) / B);
    else
        relativeError = fabs((A - B) / A);
    if (relativeError <= maxRelativeError)
        return true;
    return false;
}


#endif // INCLUDE_UTILS_H_
