//***************************************************************************
//***************************************************************************
/*
 * File: main.cpp
 * Project: euler2D
 * Description: main.cpp executes the euler2D project.  The euler2D project
 * 	illustrates the use of the simulation code to simulate
 * 	the Euler gas dynamic equations.  By default, a Sedov-like test problem
 * 	is simulated.  This project also demonstrates parallelization using MPI.
 *
 * Date: December 8, 2015
 * Author: Micah Esmond, mesmond@vt.edu
*/
//***************************************************************************
//***************************************************************************
#include <iostream>
#include <vector>
#include <ctime>
#include <chrono>
#include <cstdio>
#include "mpi.h"
#include <assert.h>
#include <string.h>

#include "dataTypes2D.h"
#include "userInputClass.h"
#include "mesmond-utils.h"
#include "structuredGeometry2D.h"
#include "fluids2D.h"


using namespace std;

void simpleEuler2D(int axialZones=100);

int main(int argc, char *argv[])
{
	//Declare
	//~ vector<PerfectGas2D> blocks;

	//Initialize the simulation based on user inputs*************************
	//~ init(argv[1], blocks);
	//~ int number_of_blocks=blocks.size();
	//~ auto wall_time_start = chrono::system_clock::now();

	int numZones[3]={100,1000,10000};
	int numMeshes=3;

	for (int i=0; i<numMeshes; ++i)
	{
		simpleEuler2D(numZones[i]);
	}
	

	return 0;
}

void simpleEuler2D(int axialZones)
{
	Point2D<double> origin(0.0,0.0);
	Point2D<double> length(1.0,0.01);
	Vector2D<int> size(axialZones,10);
	StructuredGeometry2D geom(origin, origin+length, size);
	geom.close_all_bdys();

	PerfectGas2D air(geom,1.4,1.27954e-26);


	//~ air.init_test_sedov_lowIntensity();
	air.init_test_sod_shockTube();
	air.update_boundary_values();


	Euler2D euler;
	double simTime=0.0;
	double timeStep=1.0e-17;
	double dt_dump=0.5e-6;
	double maxTime=20.0e-6; //s
	bool endSimulation=false;
	int outputCount=0;
	air.vtkOutput("output", outputCount);

	do
	{
		//~ cout << "*** SimTime=" << simTime << ", timeStep=" << timeStep << endl;
		
		for (int i=1; i<=air.getSize().get_dir1(); ++i)
		for (int j=1; j<=air.getSize().get_dir2(); ++j)
		{
			double continuity_rhs=euler.get_continuity_rhs(i,j, air);
			Vector2D<double> momentum_rhs=euler.get_momentum_rhs(i,j, air);
			double totEnergy_rhs=euler.get_totEnergy_rhs(i,j, air);

			air.continuity_write_rhs(i,j, continuity_rhs);
			air.momentum_write_rhs(i,j, momentum_rhs);
			air.totalEnergy_write_rhs(i,j, totEnergy_rhs);
		}

		for (int i=1; i<=air.getSize().get_dir1(); ++i)
		for (int j=1; j<=air.getSize().get_dir2(); ++j)
		{
			air.update_from_rhs(i,j, timeStep);
		}

		simTime+=timeStep;
		air.update_boundary_values();

		if (simTime >= outputCount*dt_dump)
			air.vtkOutput("output", outputCount);

		timeStep=min_2arg(air.get_explicit_timeStep(), 1.5*timeStep);

		if (simTime+timeStep > maxTime)
		{
			timeStep=maxTime-simTime;
			if (isEqual(0.0, timeStep))
				endSimulation=true;
		}

	} while (simTime <= maxTime && endSimulation==false);

	air.sodProblemOutput(simTime);


	cout << "Simulation Done!" << endl;
}

