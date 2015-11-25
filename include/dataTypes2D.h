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
#include <iomanip>

#include "mesmond-utils.h"
#include "simulationConstants.h"

using namespace std;

//***************************************************************************
//***************************************************************************
//***************************************************************************
//Point2D********************************************************************
/*
 * Class: Point2D
 * Purpose: Store coordinates of a particular point in (dir1, dir2).
 *
 * Description: dir1 and dir2 correspond to x and y in cartesion, or r and z
 * 	in cylindrical coordinates.  dir1 is along the  horizontal
 * 	axis, and dir2 is along the vertical axis.
 */
template <class type> class Point2D
{
	protected:
	type dir1;
	type dir2;

	public:
	explicit Point2D(type dir1=0, type dir2=0)
		: dir1(dir1), dir2(dir2) {};

	inline void write_dir1(const type& value) { dir1=value; }
	inline void write_dir2(const type& value)	{ dir2=value; }
	inline type get_dir1() const { return dir1; }
	inline type get_dir2() const { return dir2; }
};

//Define Operators***********************************************************
template <typename T> ostream& operator<<(ostream& os, Point2D<T> point)
{
	os << "(" << point.get_dir1() << "," << point.get_dir2() << ")";
	return os;
}

template <typename T> Point2D<T> operator-(const Point2D<T>& p1, const Point2D<T>& p2)
{
	return Point2D<T>(	p1.get_dir1() - p2.get_dir1(),
						p1.get_dir2() - p2.get_dir2() 	);
}

template <typename T> Point2D<T> operator+(const Point2D<T>& p1, const Point2D<T>& p2)
{
	return Point2D<T>(	p1.get_dir1() + p2.get_dir1(),
						p1.get_dir2() + p2.get_dir2() 	);
}

template <typename T> Point2D<T> operator*(const Point2D<T>& point, const T& value)
{
	return Point2D<T>(	point.get_dir1() * value,
						point.get_dir2() * value 	);
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
	explicit Vector2D(type dir1=0, type dir2=0)
		: Point2D<type>(dir1, dir2) {};

	inline type get_magnitude() const
		{ return sqrt(pow(this->get_dir1(),2.0)+pow(this->get_dir2(),2.0)); }
};


//Define Operators***********************************************************
template <typename T> Vector2D<T> operator*(const Vector2D<T>& vector, const T& value)
{
	return Vector2D<T>(	vector.get_dir1() * value,
						vector.get_dir2() * value 	);
}

template <typename T> Vector2D<T> operator*(const T& value, const Vector2D<T>& vector)
{
	return vector*value;
}

template <typename T> Vector2D<T> operator-(const Vector2D<T>& v1, const Vector2D<T>& v2)
{
	return Vector2D<T>(	v1.get_dir1() - v2.get_dir1(),
						v1.get_dir2() - v2.get_dir2() 	);
}

template <typename T> Vector2D<T> operator+(const Vector2D<T>& v1, const Vector2D<T>& v2)
{
	return Vector2D<T>(	v1.get_dir1() + v2.get_dir1(),
						v1.get_dir2() + v2.get_dir2() 	);
}


template <typename T> T dotProduct(const Vector2D<T>& vec1, const Vector2D<T>& vec2)
{
	T product0=vec1.get_dir1()*vec2.get_dir1();
	T product1=vec1.get_dir2()*vec2.get_dir2();
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
 * 	together, the result will be a Dyad2D as in normal mathematics.
 *
 * The Dyad in 2D is a 2x2 matrix.  It is output from
 * 	the dyadic product of two vectors.  Each element of a
 * 	dyad is a multiplicative arrangment of the two components of the
 * 	vectors of interest.  In this implementation, value_01 is the
 * 	product of of the zeroth component of the first vector in the
 * 	multiplication and the first component of the second vector in the
 * 	multiplication.
 */
class Dyad2D
{
	private:
	double value_00;
	double value_01;
	double value_10;
	double value_11;

	public:
	explicit Dyad2D(int fillValue=0)
	{
		value_00=(double)fillValue;
		value_01=(double)fillValue;
		value_10=(double)fillValue;
		value_11=(double)fillValue;
	}

	friend ostream& operator<<(ostream& os, const Dyad2D vector);
	friend Dyad2D operator*(const Vector2D<double>& vec0, const Vector2D<double>& vec1);
	friend Dyad2D operator*(const double& scalar, const Dyad2D& dyad);

	/*
	 * Function: get_dir1()
	 * Purpose: Returns the portion of the Dyad realated to the
	 * 	primary direction (i.e. x in (x,y), r in (r,z)).
	 *
	 * Description: If the multiplication of the vectors is assumed to
	 * 	u*v where u and v are 2D vectors, get_dir1() of the resulting
	 * 	dyad provides the "amount" the zeroth component of
	 * 	vector v carried by vector u.  Similarly, get_dir2()
	 * 	returns the "amount" amount of the first component of
	 * 	vector v carried by vector u.
	 * 	
	 */
	Vector2D<double> get_dir1() const
		{ return Vector2D<double>(value_00, value_10); }

	Vector2D<double> get_dir2() const
		{ return Vector2D<double>(value_01, value_11); }
};

ostream& operator<<(ostream& os, const Dyad2D vector);
Dyad2D operator*(const Vector2D<double>& vec0, const Vector2D<double>& vec1);
Dyad2D operator*(const double& scalar, const Dyad2D& dyad);










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
 * 	the convention that dir1 (x or r) points to the East, while
 * 	dir2 (y or z) points to the North.
 * 	The P stands for Point and is a placeholder for values at the cell
 * 	center. The directions dir1 and dir2 are included for the P placeholder
 * 	so that fluxes in these directions can be stored and used.
 */
template <class type> class StructuredLocalField2D
{
	public:
	type P;
	type P_dir1;
	type P_dir2;
	type N;
	type S;
	type E;
	type W;

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
	int numElements;
	int maxHeight;
	type *array;

	void init(Vector2D<int> size_to_init)
	{
		size=size_to_init;
		numElements=(size.get_dir1()+2)*(size.get_dir2()+2);
			//Add extra slots for bdy cells
		maxHeight=size.get_dir2()+2;
	}

	void reset()	
	{
		size=Vector2D<int>(0,0);
		numElements=0;
		maxHeight=0;
	}
	

	public:
	explicit SpacialArray2D( Vector2D<int> size=Vector2D<int>(10,10) ); //Con
	SpacialArray2D(const SpacialArray2D<type>& that); //Copy Con
	SpacialArray2D(SpacialArray2D<type>&& other); //Move Con
	SpacialArray2D<type>& operator=(SpacialArray2D<type>&& other); //Move Assign
	SpacialArray2D<type> operator=(const SpacialArray2D<type>& rhs); //Copy Assign
	~SpacialArray2D(); //Destructor

	inline type get(int i, int j) const { return array[i*maxHeight+j]; }
	inline Vector2D<int> getSize() const { return size; }
	inline int getSize_dir1() const { return size.get_dir1(); }
	inline int getSize_dir2() const { return size.get_dir2(); }
	inline int getCount_dir1() const { return size.get_dir1(); }
	inline int getCount_dir2() const { return size.get_dir2(); }
	inline void write(int i, int j, const type& value) { array[i*maxHeight+j]=value; }

	void print() const;
	void fill(const type& value);
	void set_adiabaticBdyValues();

	StructuredLocalField2D<type> getLocalField(int i, int j) const;
};

//Define Operators***********************************************************
template <typename T> SpacialArray2D<T> operator*(
	const SpacialArray2D<T>& array, const double& scalar)
{
	SpacialArray2D<T> result=array;

	for (int i=0; i<=result.getCount_dir1()+1; ++i)
	for (int j=0; j<=result.getCount_dir2()+1; ++j)
	{
		result.write(i,j, array.get(i,j)*scalar);
	}
	return result;
}

#endif  // INCLUDE_DATATYPES2D_H_
