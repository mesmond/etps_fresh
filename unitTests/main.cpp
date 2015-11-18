//***************************************************************************
//***************************************************************************
/*
 * File: main.cpp
 * Project: unitTests
 * Description: main.cpp executes the unitTests project which tests the code
 * 	under development for bugs and deficiencies.
 *
 * Date: November 17, 2015
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
#include "fluids2D.h"
#include "operators2D.h"


using namespace std;



void test_SpacialArray2D();

int main(int argc, char *argv[])
{
	cout << "************************************" << endl;
	cout << "Testing Operators2D_cyl ..." << endl;

	SpacialArray2D<double> pressure(Vector2D<int>(20,5));
	Vector2D<double> gradP;
	pressure.fill(5.0);
	pressure.print();
	

	Operators2D operate;
	gradP=operate.gradient(2,3, pressure);

	cout << "gradP=" << gradP << endl;

}

void test_PerfectGas2D()
{
	cout << "************************************" << endl;
	cout << "Testing PerfectGas2D ..." << endl;
	PerfectGas2D air(Vector2D<int>(15,30));

	double air_temperature=20;	//deg C
	air.fill_temperature(air_temperature+273.15);	//K
	air.fill_pressure(101325.0);					//Pa
	air.fill_velocity(Vector2D<double>(0.0,0.0));	//m/s
	air.init_fromBasicProps();


	
	int soundSpeed=(int)air.getSoundSpeed(2,3);
	cout << "soundSpeed=" << soundSpeed << endl;
	cout << "thermalConductivity=" << air.getThermalConductivity(2,3) << endl;
	
	assert(soundSpeed == 343); // m/s
	cout << "Done***************************************" << endl;
}


void test_SpacialArray2D()
{
	cout << "************************************" << endl;
	cout << "Testing SpacialArray2D ..." << endl;

	SpacialArray2D< Vector2D<double> > velocity(Vector2D<int>(10,5));
	Vector2D<double> localVel(4.5,3.4);

	velocity.fill(Vector2D<double>(2.0,3.0));

	velocity.write(2,2, localVel);
	velocity.print();

	StructuredLocalField2D<Vector2D<double> > localField=velocity.getLocalField(2,1);
	localField.print();	

	assert( isEqual(velocity.get(2,2).get_dir0(), localField.N.get_dir0()) );
	assert( isEqual(velocity.get(2,0).get_dir0(), localField.S.get_dir0()) );
	assert( isEqual(velocity.get(3,1).get_dir0(), localField.E.get_dir0()) );
	assert( isEqual(velocity.get(1,1).get_dir0(), localField.W.get_dir0()) );
}

void test_Point2D_Vector2D()
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

	cout << "------------------------------------" << endl;
	cout << "Vector Magnitude--------------------" << endl;
	vec=Vector2D<double>(3.0,4.0);
	cout << "vec       =" << vec << endl;
	cout << "vec (mag) =" << vec.get_magnitude() << endl;
	assert( isEqual(vec.get_magnitude(), 5.0) );

	Vector2D<int> vecInt(3,4);
	cout << "vecInt      =" << vecInt << endl;
	cout << "vecInt (mag)=" << vecInt.get_magnitude() << endl;
	assert( vecInt.get_magnitude() == 5 );

	cout << "------------------------------------" << endl;
	cout << "Vector DotProduct-------------------" << endl;
	vec=Vector2D<double>(3.0,4.0);
	vec2=Vector2D<double>(3.0,4.0);
	double dotProd=dotProduct(vec, vec2);

	cout << "vec        =" << vec << endl;
	cout << "vec2       =" << vec2 << endl;
	cout << "dotProd    =" << dotProd << endl;

	assert( isEqual(dotProd, 25.0) );
}


void test_Dyad2D()
{
	cout << "************************************" << endl;
	cout << "Testing Dyad2D ..." << endl;

	Vector2D<double> vec(4,6);
	Vector2D<double> vec2(3,5);

	Dyad2D dyad;

	cout << "vec =" << vec << endl;
	cout << "vec2=" << vec2 << endl;
	
	cout << "dyad (init)=" << dyad << endl;
	dyad=vec*vec2;
	cout << "dyad=vec*vec2 : dyad=" << dyad << endl;
}



void test_StructuredGeometry2D()
{
	cout << "************************************" << endl;
	cout << "Testing StructuredGeometry2D ..." << endl;

	Point2D<double> origin(0.0,0.0);
	Point2D<double> extent(1.0,3.5);
	Vector2D<int> size(34,34);
	
	StructuredGeometry2D block(origin,extent,size);

	origin=Point2D<double>(0.0,3.5);
	extent=Point2D<double>(1.0,4.5);
	size=Vector2D<int>(34,15);
	
	StructuredGeometry2D block_north(origin,extent,size);

	origin=Point2D<double>(1.0,3.5);
	Point2D<double> length(1.0,3.25);
	size=Vector2D<int>(34,15);
	
	StructuredGeometry2D block_east(origin,origin+length,size);

	block.link_north(block_north);
	block_north.link_south(block);

	block_north.link_east(block_east);
	block_east.link_west(block_north);
}
