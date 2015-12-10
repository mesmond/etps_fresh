//***************************************************************************
//***************************************************************************
/*
 * File: mesmond-utils.cpp
 * Description: This file contains definitions of
 * 	several useful one-off functions. It includes the folowing:
 * 		min_2arg (double)
 * 		max_array (double)
 * 		isEqual (double)
*/
//***************************************************************************
//***************************************************************************
#include "mesmond-utils.h"

double min_2arg(double arg1, double arg2)
{
	double min=arg1;
	if (arg2<min) min=arg2;
	return min;
}

double max_2arg(double arg1, double arg2)
{
	double max=arg1;
	if (arg2>max) max=arg2;
	return max;
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


int getLength(const char* word)
{
	int length=0;
	char thisChar=word[length];
	while (thisChar != '\0')
	{
		length++;
		thisChar=word[length];
	}

	return length;
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
bool isEqual(double A, double B, double maxRelativeError)
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
