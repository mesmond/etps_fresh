//***************************************************************************
//***************************************************************************
/*
 * File: Utils.h
 * Description: This file contains declarations and definitions of
 * 	several useful functions.
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


void outputToLogFile(
	string fileName,
	double simTime,
	double timeStep,
	int timeStepCount,
	double wallTime)
{
	ofstream output;

	output.open(fileName.c_str(), ios::app);

	if (output.is_open())
	{
		output << "************************************" << endl;
		output << "wallTime(s)  =" << wallTime << endl;
		output << "simTime(s)   =" << simTime << endl;
		output << "timeStep(s)  =" << timeStep << endl;
		output << "timeStepCount=" << timeStepCount << endl;
	}
	else
	{
		cout << "ERROR: Unable to open Log file!!!" << endl;
		exit(1);
	}

}


#endif // INCLUDE_UTILS_H_
