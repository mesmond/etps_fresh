//***************************************************************************
//***************************************************************************
/*
 * File: dataTypes2D.h
 * Description: This file contains class declarations for 2D data types.
 * 	The datatypes included are:
 * 		Point2D
 * 		Vector2D
 * 		Dyad2D
 * 		StructuredLocalField2D
 * 		StructuredGeometry2D
 * 		SpacialArray2D
*/
//***************************************************************************
//***************************************************************************
#ifndef INCLUDE_DATATYPES2D_H_
#define INCLUDE_DATATYPES2D_H_

#include <iostream>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "mesmond-utils.h"

using namespace std;

//***************************************************************************
//***************************************************************************
//***************************************************************************
//Point2D********************************************************************
/*
 * Class: Point2D
 * Purpose: Store coordinates of a particular point in (dir0, dir1).
 *
 * Description: dir0 and dir1 correspond to x and y in cartesion, or r and z
 * 	in cylindrical coordinates.  dir0 is along the  horizontal
 * 	axis, and dir1 is along the vertical axis.
 */
template <class type> class Point2D
{
	protected:
	type dir0;
	type dir1;

	public:
	explicit Point2D(type dir0=0, type dir1=0)
		: dir0(dir0), dir1(dir1) {};

	inline void write_dir0(const type& value) { dir0=value; }
	inline void write_dir1(const type& value)	{ dir1=value; }
	inline type get_dir0() const { return dir0; }
	inline type get_dir1() const { return dir1; }
};

//Define Operators***********************************************************
template <typename T> ostream& operator<<(ostream& os, Point2D<T> point)
{
	os << "(" << point.get_dir0() << "," << point.get_dir1() << ")";
	return os;
}

template <typename T> Point2D<T> operator-(const Point2D<T>& p1, const Point2D<T>& p2)
{
	Point2D<T> result;

	result.write_dir0(p1.get_dir0() - p2.get_dir0());
	result.write_dir1(p1.get_dir1() - p2.get_dir1());

	return result;
}

template <typename T> Point2D<T> operator+(const Point2D<T>& p1, const Point2D<T>& p2)
{
	Point2D<T> result;

	result.write_dir0(p1.get_dir0() + p2.get_dir0());
	result.write_dir1(p1.get_dir1() + p2.get_dir1());

	return result;
}

template <typename T> Point2D<T> operator*(const Point2D<T>& point, const T& value)
{
	Point2D<T> result;
	result.write_dir0(point.get_dir0() * value);
	result.write_dir1(point.get_dir1() * value);
	return result;
}

template <typename T> Point2D<T> operator*(const T& value, const Point2D<T>& point)
{
	return point*value;
}








//***************************************************************************
//***************************************************************************
//***************************************************************************
//Vector2D*******************************************************************
/*
 * Class: Vector2D
 * Purpose: Act as a mathematical vector in 2D.
 *
 * Description: This class simply stores and processes a mathematical
 * 	vector.  This is not to be confused with a std::vector.	Vector2D
 * 	utilizes the functionality of a Point2D type while adding
 * 	a get_magnitude() function.
 */
template <class type> class Vector2D : public Point2D<type>
{	

	public:
	explicit Vector2D(type dir0=0, type dir1=0)
		: Point2D<type>(dir0, dir1) {};

	inline type get_magnitude() const
	{
		return sqrt(
			 pow(this->get_dir0(),2.0)
			+pow(this->get_dir1(),2.0));
	}
};


//Define Operators***********************************************************
template <typename T> T dotProduct(const Vector2D<T>& vec1, const Vector2D<T>& vec2)
{
	T product0=vec1.get_dir0()*vec2.get_dir0();
	T product1=vec1.get_dir1()*vec2.get_dir1();

	return product0+product1;
}










//***************************************************************************
//***************************************************************************
//***************************************************************************
//Dyad2D*********************************************************************
/*
 * Class: Dyad2D
 * Purpose: Act as a mathematical Dyad in 2D.
 *
 * Description: This class stores and processes Dyad objects in 2D.  Only a
 * 	double type is allowed for its members.  When two vectors are multiplied
 * 	together, the result will be this Dyad2D as in normal mathematics.
 */
class Dyad2D
{
	private:
	double value_00;
	double value_01;
	double value_10;
	double value_11;

	public:
	Dyad2D()
	{
		value_00=0.0;
		value_01=0.0;
		value_10=0.0;
		value_11=0.0;
	}

	friend ostream& operator<<(ostream& os, const Dyad2D vector);
	friend Dyad2D operator*(const Vector2D<double>& vec0, const Vector2D<double>& vec1);
};

ostream& operator<<(ostream& os, const Dyad2D vector);
Dyad2D operator*(const Vector2D<double>& vec0, const Vector2D<double>& vec1);











//***************************************************************************
//***************************************************************************
//***************************************************************************
//StructuredLocalField2D*****************************************************
/*
 * Class: StructuredLocalField2D
 * Purpose: Provide a basic object for acquiring and sending simulation
 * 	information relevant to a particular cell.
 *
 * Description: North, South, East, and West are determined according to
 * 	the convention that dir0 (x or r) points to the East, while
 * 	dir1 (y or z) points to the North.
 * 	The P stands for Point and is a placeholder for values at the cell
 * 	center. The directions dir0 and dir1 are included for the P placeholder
 * 	so that fluxes in these directions can be stored and used.
 */
template <class type> class StructuredLocalField2D
{
	public:
	type P;
	type P_dir0;
	type P_dir1;
	type N;
	type S;
	type E;
	type W;

	StructuredLocalField2D()
	{
		P=(type)0;
		P_dir0=(type)0;
		P_dir1=(type)0;
		N=(type)0;
		S=(type)0;
		E=(type)0;
		W=(type)0;
	}

	void print() const;
};












//***************************************************************************
//***************************************************************************
//***************************************************************************
//SpacialArray2D*************************************************************
/*
 * Class: SpacialArray2D
 * Purpose: Provide a template class to store and process simulated
 * 	spacial distriutions (e.g. pressure, velocity, etc.).
 *
 * Description: This template class allows the user to define
 * 	a spacial array of different types such as a vector (Vector2D)
 * 	or a scalar (double).
 * 	This class provides an extension of the basic 2D data types to include
 * 	spacial distributions.
 */
template <class type> class SpacialArray2D
{
	private:
	Vector2D<int> size;
	type **array;

	public:
	explicit SpacialArray2D( Vector2D<int> size=Vector2D<int>(10,10) );
	SpacialArray2D(const SpacialArray2D<type>& that);
	~SpacialArray2D();

	inline type get(int i, int j) const { return array[i][j]; }
	inline int getSize_dir0() const { return size.get_dir0(); }
	inline int getSize_dir1() const { return size.get_dir1(); }
	inline int getCount_dir0() const { return size.get_dir0(); }
	inline int getCount_dir1() const { return size.get_dir1(); }
	inline void write(int i, int j, const type& value) { array[i][j]=value; }

	void print() const;

	void fill(const type& value);

	StructuredLocalField2D<type> getLocalField(int i, int j);
};














//***************************************************************************
//***************************************************************************
//***************************************************************************
//StructuredGeometry2D*******************************************************
/*
 * class: StructuredGeometry2D
 * Purpose: Hold and process geometric parameters for a structured geometry.
 *
 * Description: StructuredGeometry2D holds information on the spacial placement
 * 	of a particular geometry. It also holds information
 * 	regarding the mesh spacing and mesh density.
 * 	Global Coordinates of structured cells
 * 	are specified using a cell-centered scheme.
 * 	This struct is designed for an orthogonal grid.
 */
class StructuredGeometry2D
{
	private:
	Point2D<double> origin;
	Point2D<double> extent;
	Vector2D<double> length;

	Vector2D<int> numZones;

	struct MeshSpacing
	{
		double *dir0;
		double *dir1;
	} zoneDelta, globalCoord;

	double minMeshSpacing;

	public:
	StructuredGeometry2D(const StructuredGeometry2D& that);
	explicit StructuredGeometry2D(
		Point2D<double> origin=Point2D<double>(0.0,0.0),
		Point2D<double> extent=Point2D<double>(1.0,1.0),
		Vector2D<int> size=Vector2D<int>(10,10),
		double refineMesh_dir0=1.0, double refineMesh_dir1=1.0);
	~StructuredGeometry2D();

	//Get Information********************************************************
	inline double get_zoneDelta_north() const
		{ return zoneDelta.dir1[numZones.get_dir1()]; }

	inline double get_zoneDelta_south() const
		{ return zoneDelta.dir1[1]; }

	inline double get_zoneDelta_east() const
		{ return zoneDelta.dir0[numZones.get_dir0()]; }

	inline double get_zoneDelta_west() const { return zoneDelta.dir0[1]; }

	inline int getCount_dir0() const { return numZones.get_dir0(); }
	inline int getCount_dir1() const { return numZones.get_dir1(); }
	inline Point2D<double> get_origin() const { return origin; }
	inline Point2D<double> get_extent() const { return extent; }


	//Mutation Handling******************************************************
	void link_north(const StructuredGeometry2D& toLink);
	void link_south(const StructuredGeometry2D& toLink);
	void link_east(const StructuredGeometry2D& toLink);
	void link_west(const StructuredGeometry2D& toLink);


	private:
	//Set Information********************************************************
	inline void setBdy_zoneDelta_north(double value)
	{
		int index=numZones.get_dir1()+1;
		zoneDelta.dir1[index] = value;
		globalCoord.dir1[index]=extent.get_dir1()+0.5*zoneDelta.dir1[index];
	}

	inline void setBdy_zoneDelta_south(double value)
	{
		int index=0;
		zoneDelta.dir1[index] = value;
		globalCoord.dir1[index]=origin.get_dir1()-0.5*zoneDelta.dir1[index];
	}

	inline void setBdy_zoneDelta_east(double value)
	{
		int index=numZones.get_dir0()+1;
		zoneDelta.dir0[index] = value;
		globalCoord.dir0[index]=extent.get_dir0()+0.5*zoneDelta.dir0[index];
	}

	inline void setBdy_zoneDelta_west(double value)
	{
		int index=0;
		zoneDelta.dir0[index] = value;
		globalCoord.dir0[index]=origin.get_dir0()-0.5*zoneDelta.dir0[index];
	}


};



















//~ 
//~ 
//~ 
//~ 
//~ 
//~ /*
 //~ * Class: spacialArray2D
 //~ * Purpose: Provide a basic object for spacial arrays.
 //~ */
//~ class spacialArray2D
//~ {
	//~ protected:
	//~ struct arraySize
	//~ {
		//~ int dir0;
		//~ int dir1;
	//~ } numZones;
//~ 
	//~ double **array;
//~ 
	//~ public:
	//~ spacialArray2D(const spacialArray2D& that);
	//~ spacialArray2D(int size_dir0=10, int size_dir1=10, double valueInit=0.0);
	//~ spacialArray2D& operator=(const spacialArray2D& that);
//~ 
	//~ 
	//~ ~spacialArray2D();
//~ 
	//~ double get(int i, int j) const;
	//~ void write(int i, int j, double value);
//~ 
	//~ int getSize_dir0() const;
	//~ int getSize_dir1() const;
//~ 
	//~ cellDirections getLocalField(int i, int j) const;
//~ };
//~ 
//~ /*
 //~ * class: scalarArray2D
 //~ * Purpose: Store and process 2D scalar arrays.
 //~ */
//~ class scalarArray2D : public spacialArray2D
//~ {
	//~ public:
	//~ scalarArray2D(const scalarArray2D& that) : spacialArray2D(that){}
	//~ scalarArray2D(int size_dir0=10, int size_dir1=10)
		//~ : spacialArray2D(size_dir0, size_dir1) {}
//~ };
//~ 
//~ /*
 //~ * class: vectorArray2D
 //~ * Purpose: Store and process 2D vector arrays.
 //~ */
//~ class vectorArray2D
//~ {
	//~ private:
	//~ spacialArray2D array_dir0;
	//~ spacialArray2D array_dir1;
//~ 
	//~ public:
	//~ vectorArray2D(const vectorArray2D& that)
		//~ : 	array_dir0(that.array_dir0),
			//~ array_dir1(that.array_dir1) {}
	//~ vectorArray2D(int size_dir0=10, int size_dir1=10)
		//~ :	array_dir0(size_dir0, size_dir1),
			//~ array_dir1(size_dir0, size_dir1) {}
//~ 
	//~ double getMagnitude(int i, int j);
	//~ void writeComponent_dir0(int i, int j, double value);
	//~ void writeComponent_dir1(int i, int j, double value);
//~ 
	//~ double component_dir0(int i, int j);
	//~ double component_dir1(int i, int j);
//~ 
	//~ int getSize_dir0();
	//~ int getSize_dir1();
//~ 
	//~ spacialArray2D dir0() {return array_dir0;}
//~ 
	//~ spacialArray2D dir1() {return array_dir1;}
//~ 
	//~ cellDirections product(cellDirections scalarField, int i, int j) const;
//~ 
	//~ cellDirections localField_dir0(int i, int j)
	//~ {
		//~ return array_dir0.getLocalField(i,j);
	//~ }
	//~ cellDirections localField_dir1(int i, int j)
	//~ {
		//~ return array_dir1.getLocalField(i,j);
	//~ }
//~ };


#endif  // INCLUDE_DATATYPES2D_H_
