//***************************************************************************
//***************************************************************************
/*
 * File: main.cpp
 * Project: euler2D
 * Description: main.cpp executes the euler2D simulation code. The simulation
 * 	domain is divided into several regions.  The data for each region is
 * 	stored in a MeshBlock object and processed by an individual MPI process.
 * 	MeshBlock objects are stored in a single vector.
 *
 * 	An init() function is used to initialize the simulation variables.
 * 	The governing equations are used to advance the simulation in
 * 	euler2D_timeIntegration().
 *
 * 	The simulation wall time is tracked and information is written to a log
 * 	file periodically.  Output files are announced in the terminal console.
 *
 * Date: October 7, 2015
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


#include "dataTypes2D.h"
#include "userInputClass.h"
#include "mesmond-utils.h"

using namespace std;

int main(int argc, char *argv[])
{
	cout << "************************************" << endl;
	cout << "Testing Point2D and Vector2D..." << endl;

	Point2D<double> point(3.0, 4.6);
	Vector2D<double> vec(4.5,3.6);
	double scale=3.0;
	
	cout << "point=" << point << endl;
	cout << "vec  =" << vec << endl;
	
	cout << "point-vec=" << point-vec << endl;
	cout << "point+vec=" << point+vec << endl;

	cout << "------------------------------------" << endl;
	Point2D<double> point2(7.8,76.9);
	point2=vec*scale;
	cout << "scale=" << scale << endl;
	cout << "vec*scale=" << vec	*scale << endl;
	cout << "point*scale=" << point*scale << endl;
	assert( isEqual(point2.get_dir1(), vec.get_dir1()*scale) );
	assert( isEqual(point2.get_dir0(), vec.get_dir0()*scale) );


	cout << "------------------------------------" << endl;
	Vector2D<double> vec2(3.5,2.5);
	vec=vec2;
	cout << "vec2 =" << vec2 << endl;
	cout << "vec=vec2:" << endl;
	cout << "vec =" << vec << endl;
	cout << "vec2=" << vec2 << endl;
	assert( isEqual(vec.get_dir1(), vec2.get_dir1()) );
	assert( isEqual(vec.get_dir0(), vec2.get_dir0()) );

	
	
	//~ StructuredGeometry2D block(0.0,0.0,1.0,1.0,10,10);
	//~ StructuredGeometry2D block_north(0.0,1.0,1.0,3.0,10,5);
//~ 
	//~ block.link_north(block_north);
	//~ block_north.link_south(block);
//~ 
	//~ cout << "block_north:" << endl;
	//~ for (int i=block_north.numZones.dir1+1; i>=0; --i)
	//~ {
		//~ ////~ cout << block_north.zoneDelta.dir1[i] << endl;
		//~ cout << block_north.globalCoord.dir1[i] << endl;
	//~ }
//~ 
	//~ 
	//~ cout << "block:" << endl;
	//~ for (int i=block.numZones.dir1+1; i>=0; --i)
	//~ {
		//~ ////~ cout << block.zoneDelta.dir1[i] << endl;
		//~ cout << block.globalCoord.dir1[i] << endl;
	//~ }



}

void testVector2D()
{
	Vector2D<double> vec(3,4);
	Vector2D<double> vec2(4,5);

	Vector2D<double> velocity(6,2);

	double scalarVal=3.0;

	cout << "vec =" << vec << endl;
	cout << "vec2=" << vec2 << endl;
	cout << "vec dot vec2=" << dotProduct(vec, vec2) << endl;
	cout << "vec*vec2=" << vec*vec2 << endl;
	cout << "vec2 - vec=" << vec2 - vec << endl;
	cout << "vec - vec2=" << vec - vec2 << endl;

	cout << "scalarVal*vec=" << scalarVal*vec << endl;
	cout << "vec*scalarVal=" << vec*scalarVal << endl;


	cout << "velocity=" << velocity << endl;
	cout << "Vector2D<double>(4,5) : " << Vector2D<double>(4,5) << endl;



	cout << "*****************" << endl;
	Vector2D<double> u(4,3);
	Vector2D<double> v(5,6);
	Dyad2D dyad=u*v;

	cout << "u   =" << u << endl;
	cout << "v   =" << v << endl;
	cout << "dyad=" << dyad << endl;
}
