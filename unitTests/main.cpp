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
#include <cstdio>
#include <assert.h>
#include <string.h>

#include "dataTypes2D.h"
#include "userInputClass.h"
#include "mesmond-utils.h"
#include "structuredGeometry2D.h"
#include "fluids2D.h"

using namespace std;


void test_Point2D_Vector2D();
void test_Dyad2D();
void test_SpacialArray2D();
void test_StructuredGeometry2D();
void test_PerfectGas2D();

int main(int argc, char *argv[])
{
	test_Point2D_Vector2D();
	test_Dyad2D();
	test_SpacialArray2D();
	test_StructuredGeometry2D();
	test_PerfectGas2D();
	return 0;
}


void test_PerfectGas2D()
{
	cout << "************************************" << endl;
	cout << "Testing PerfectGas2D ..." << endl;

	Point2D<double> origin(0.0,0.0);
	Point2D<double> extent(1.0,10.0);
	Vector2D<int> size(50,50);
	StructuredGeometry2D geom(origin, extent, size);

	//Test for a memory Leak.
	int i=0;
	while (i < 10)
	{
		PerfectGas2D air(geom);

		double air_temperature=20.0;	//deg C
		air.fill_temperature(air_temperature+273.15);	//K
		air.fill_pressure(101325.0);					//Pa
		air.fill_velocity(Vector2D<double>(0.0,0.0));	//m/s
		i++;
	}

	PerfectGas2D air(geom);

	double air_temperature=20.0;	//deg C
	air.fill_temperature(air_temperature+273.15);	//K
	air.fill_pressure(101325.0);					//Pa
	air.fill_velocity(Vector2D<double>(0.0,0.0));	//m/s
	air.init_from_temperature_pressure_velocity();
	//~ air.print_geometry();

	int soundSpeed=(int)air.calcSoundSpeed(2,3);
	cout << "soundSpeed=" << soundSpeed << endl;
	cout << "thermalConductivity=" << air.calcThermalConductivity(2,3) << endl;
	assert(soundSpeed == 343); // m/s



	//Cylindrical vs. Rectangular Geometry***********************************
	i=2;
	int j=3;
	origin=Point2D<double>(0.0,0.0);
	extent=Point2D<double>(1.0,1.0);
	size=Vector2D<int>(10,10);

	//Cylindrical Geometry***************************************************
	StructuredCylGeometry2D geom_cyl(origin, extent, size);
	PerfectGas2D air_cyl(geom_cyl);
	cout << "air_cyl.getVolume(" << i << "," << j << ")="
		<< air_cyl.getVolume(i,j) << endl;

	//Rectangular Geometry***************************************************
	StructuredGeometry2D geom_rec(origin, extent, size);
	PerfectGas2D air_rec(geom_rec);
	cout << "air_rec.getVolume(" << i << "," << j << ")="
		<< air_rec.getVolume(i,j) << endl;


	
	assert(sizeof(geom_cyl) == sizeof(geom_rec));
	assert( !isEqual(air_cyl.getVolume(i,j) , air_rec.getVolume(i,j)) );



	//Copy Constructor Test**************************************************
	cout << "Copy Constructor..." << endl;
	PerfectGas2D air_test(air_cyl);

	assert( isEqual(air_cyl.getVolume(i,j) , air_test.getVolume(i,j)) );
	assert(sizeof(air_cyl) == sizeof(air_test));

	PerfectGas2D air_test_rec(air_rec);

	assert( isEqual(air_test_rec.getVolume(i,j) , air_rec.getVolume(i,j)) );
	assert(sizeof(air_test_rec) == sizeof(air_rec));
	
	cout << "*******************************************" << endl;
	cout << "Done***************************************" << endl;
}


void test_SpacialArray2D()
{
	cout << "************************************" << endl;
	cout << "Testing SpacialArray2D ..." << endl;

	Vector2D<int> size(10,5);

	cout << "Creating velocity..." << endl;
	SpacialArray2D< Vector2D<double> > velocity(size);
	Vector2D<double> localVel(4.5,3.4);

	velocity.fill(Vector2D<double>(2.0,3.0));

	velocity.write(2,2, localVel);
	velocity.print();

	StructuredLocalField2D<Vector2D<double> > localField=velocity.getLocalField(2,1);
	localField.print();	

	assert( isEqual(velocity.get(2,2).get_dir1(), localField.N.get_dir1()) );
	assert( isEqual(velocity.get(2,0).get_dir1(), localField.S.get_dir1()) );
	assert( isEqual(velocity.get(3,1).get_dir1(), localField.E.get_dir1()) );
	assert( isEqual(velocity.get(1,1).get_dir1(), localField.W.get_dir1()) );


	//Test Assignment********************************************************
	SpacialArray2D<Vector2D<double> > testVec(size);
	double scalar=3.0;



	testVec=velocity*scalar;

	//~ cout << "print (velocity*scalar)..." << endl;
	//~ (velocity*scalar).print();
//~ 
	//~ cout << "print testVec..." << endl;
	//~ testVec.print();

	for (int i=0; i<=testVec.getCount_dir1()+1; ++i)
		for (int j=0; j<=testVec.getCount_dir2()+1; ++j)
		{
			//~ double value1=testVec.get(i,j).get_dir1();
			//~ double value2=velocity.get(i,j).get_dir1()*scalar;
			
			//~ cout << "Checking stuff, i=" << i << ", j=" << j;
			//~ cout << ", value1=" << value1 << ", value2=" << value2 << endl;
			assert( isEqual( testVec.get(i,j).get_dir1() ,
				velocity.get(i,j).get_dir1()*scalar) );
		}



	
	


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
	assert( isEqual(point2.get_dir2(), vec.get_dir2()*scale) );
	assert( isEqual(point2.get_dir1(), vec.get_dir1()*scale) );


	cout << "------------------------------------" << endl;
	Vector2D<double> vec2(3.5,2.5);
	vec=vec2;
	cout << "vec2 =" << vec2 << endl;
	cout << "vec=vec2:" << endl;
	cout << "vec =" << vec << endl;
	cout << "vec2=" << vec2 << endl;
	assert( isEqual(vec.get_dir2(), vec2.get_dir2()) );
	assert( isEqual(vec.get_dir1(), vec2.get_dir1()) );

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
	double scalar=5.0;

	cout << "vec    =" << vec << endl;
	cout << "vec2   =" << vec2 << endl;
	cout << "scalar =" << scalar << endl;
	
	cout << "dyad (init)=" << dyad << endl;
	
	dyad=vec*vec2;
	cout << "dyad=vec*vec2 : dyad=" << dyad << endl;

	dyad=scalar*dyad;
	cout << "dyad=scalar*dyad : dyad=" << dyad << endl;

	cout << "dyad.get_dir1()=" << dyad.get_dir1() << endl;
	cout << "dyad.get_dir2()=" << dyad.get_dir2() << endl;


	assert( isEqual(dyad.get_dir2().get_dir1(),
		vec.get_dir1()*vec2.get_dir2()*scalar));
}



void test_StructuredGeometry2D()
{
	cout << "************************************" << endl;
	cout << "Testing StructuredGeometry2D ..." << endl;

	//Testing Block Linking Functions****************************************
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



	//Test Rectangular vs. Cylindrical***************************************
	int i=2, j=3;
	origin=Point2D<double>(0.0,3.5);
	extent=Point2D<double>(1.0,4.5);
	size=Vector2D<int>(34,15);
	
	StructuredGeometry2D rec(origin, extent, size);
	StructuredCylGeometry2D cyl(origin, extent, size);

	assert( sizeof(rec) == sizeof(cyl) );
	assert( !isEqual(cyl.getVolume(i,j), rec.getVolume(i,j) ) );
	
	//~ StructuredLocalField2D<double> localField_rec=rec.getCellAreas(i,j);
	//~ StructuredLocalField2D<double> localField_cyl=cyl.getCellAreas(i,j);


	//Test Differential Operators********************************************
	size=Vector2D<int>(10,10); //Do not change
	SpacialArray2D<double> pressure(size);
	SpacialArray2D<double> divergence(size);
	SpacialArray2D<Vector2D<double> > gradP(size);
	StructuredCylGeometry2D operate(
		Point2D<double>(0.0,0.0),
		Point2D<double>(10.0,5.0),
		size);

	//Write Pressure
	for (i=1; i<=pressure.getCount_dir1(); ++i)
	for (j=1; j<=pressure.getCount_dir2(); ++j)
	{
		pressure.write(i,j, i*j);
	}

	pressure.set_NeumannBdyValues_all();

	//Compute Gradient
	for (i=1; i<=pressure.getCount_dir1(); ++i)
	for (j=1; j<=pressure.getCount_dir2(); ++j)
	{
		gradP.write( i,j, operate.gradient(i,j, pressure) );
	}

	for (i=1; i<=pressure.getCount_dir1(); ++i)
	for (j=1; j<=pressure.getCount_dir2(); ++j)
	{
		divergence.write( i,j, operate.divergence(i,j, pressure, gradP) );
	}

	operate.vtkOutput(pressure, gradP, divergence);
	
	i=3;
	j=3;
	Vector2D<double> grad=operate.gradient(i,j, pressure);
	double div1=operate.divergence(i,j, pressure, gradP);
	Vector2D<double> div2=operate.divergence(i,j, pressure, gradP, gradP*3.0);

	cout << "grad                      =" << grad << endl;
	cout << "div1(scalar,vector)       =" << div1 << endl;
	cout << "div2(scalar,vector,vector)=" << div2 << endl;

	assert( isEqual(div1, 55.8) );
	assert( isEqual(div2.get_dir1(), 826.2) );
	assert( isEqual(div2.get_dir2(), 1177.2) );
	
	assert( isEqual(grad.get_dir1(), 3.0) );
	assert( isEqual(grad.get_dir2(), 6.0) );
	


	cout << "*******************************************" << endl;
	cout << "Done***************************************" << endl;
}
