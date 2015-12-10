//***************************************************************************
//***************************************************************************
/*
 * File: dataTypes2D.cpp
 * Description: This file contains class definitions for 2D array data types.
*/
//***************************************************************************
//***************************************************************************
#include "dataTypes2D.h"

//***************************************************************************
//***************************************************************************
//***************************************************************************
//Point2D********************************************************************




//***************************************************************************
//***************************************************************************
//***************************************************************************
//Vector2D*******************************************************************

//Define Operators***********************************************************
bool operator==(const Vector2D<int>& vec1, const Vector2D<int>& vec2)
{
	if ( vec1.get_dir1() != vec2.get_dir1() )
		return false;
	else if ( vec1.get_dir2() != vec2.get_dir2() )
		return false;

	return true;
}



//***************************************************************************
//***************************************************************************
//***************************************************************************
//Dyad2D*********************************************************************

//Define Operators***********************************************************
ostream& operator<<(ostream& os, const Dyad2D dyad)
{
	os << "[(" << dyad.value_00 << "," << dyad.value_01 << "),("
		<< dyad.value_10 << "," << dyad.value_11 << ")]";

	return os;
}

Dyad2D operator*(const Vector2D<double>& vec0, const Vector2D<double>& vec1)
{
	Dyad2D dyad;

	dyad.value_00=vec0.get_dir1()*vec1.get_dir1();
	dyad.value_01=vec0.get_dir1()*vec1.get_dir2();
	dyad.value_10=vec0.get_dir2()*vec1.get_dir1();
	dyad.value_11=vec0.get_dir2()*vec1.get_dir2();

	return dyad;
}

Dyad2D operator*(const double& scalar, const Dyad2D& dyad)
{
	Dyad2D result;

	result.value_00=dyad.value_00*scalar;
	result.value_01=dyad.value_01*scalar;
	result.value_10=dyad.value_10*scalar;
	result.value_11=dyad.value_11*scalar;

	return result;
}

//***************************************************************************
//***************************************************************************
//***************************************************************************
//StructuredLocalField2D*****************************************************
template <typename type> void StructuredLocalField2D<type>::print() const
{
	cout << "StructuredLocalField2D:" << endl;
	cout << "\t P      =" << P << endl;
	cout << "\t P(dir1)=" << P_dir1 << endl;
	cout << "\t P(dir2)=" << P_dir2 << endl;
	cout << "\t N      =" << N << endl;
	cout << "\t S      =" << S << endl;
	cout << "\t E      =" << E << endl;
	cout << "\t W      =" << W << endl;
}

template class StructuredLocalField2D<Vector2D<double> >;
template class StructuredLocalField2D<Dyad2D >;
template class StructuredLocalField2D<double>;
template class StructuredLocalField2D<int>;

//***************************************************************************
//***************************************************************************
//***************************************************************************
//SpacialArray2D*************************************************************
template <typename type> SpacialArray2D<type>::SpacialArray2D( Vector2D<int> size )
{
	//Constructor
	//~ cout << "In the constructor..." << endl;
	init(size);

	array = new type [numElements];
	for (int i=0; i<numElements; ++i) 
		array[i]=(type)0;

}

template <typename type> SpacialArray2D<type>::SpacialArray2D(const SpacialArray2D<type>& that)
{
	//Copy Constructor
	//~ cout << "In the copy constructor" << endl;
	init(that.size);

	array = new type [numElements];
	for (int i=0; i<numElements; ++i)
		array[i]=that.array[i];
}

template <typename type> SpacialArray2D<type>::SpacialArray2D(SpacialArray2D<type>&& other)
	:	size(0),
		array(nullptr)
{
	//Move Constructor
	cout << "In the Move Constructor..." << endl;

	init(other.size);
	array=other.array;

	other.array=nullptr;
	other.reset();
}

template <typename type> SpacialArray2D<type>&
	SpacialArray2D<type>::operator=(SpacialArray2D<type>&& other)
{
	//Move Assignment Operator
	cout << "In the Move Assignment Operator..." << endl;

	if (this != &other)
	{
		//Delete Memory of the present object.
		delete[] array;
		array=nullptr;
		reset();

		//Transfer Data
		init(other.size);
		array=other.array;

		//Reset the other objects pointers to
		//	keep the memory from being freed
		//	by the destructor.
		other.array=nullptr;
		other.reset();
	}
	return *this;
}


template <typename type> SpacialArray2D<type>
	SpacialArray2D<type>::operator=(const SpacialArray2D<type>& rhs)
{
	//Assignment Operator
	cout << "In assignment operator.. " << endl;
	if (this != &rhs)
	{
		assert(size==rhs.size);

		for (int i=0; i<numElements; ++i)
			array[i]=rhs.array[i];
	}
	return *this;
}

template <typename type> SpacialArray2D<type>::~SpacialArray2D()
{
	//Destructor
	//~ cout << "In the Destructor!*********************************" << endl;

	if (array != nullptr)
	{
		delete[] array;
		array=nullptr;
		reset();
	}
	else
	{
		cout << "null ptr detected before free in SpacialArray2D destuctor." << endl;
	}
}


template <typename type> type SpacialArray2D<type>::getSum()
{
	type summation=(type)0;
	for (int i=1; i<=getSize().get_dir1(); ++i)
	for (int j=1; j<=getSize().get_dir2(); ++j)
	{
		summation=summation+get(i,j);
	}

	return summation;
}

template <typename type> void SpacialArray2D<type>::print() const
{
	for (int j=size.get_dir2()+1; j>=0; --j)
	{
		for (int i=0; i<=size.get_dir1()+1; ++i)
		{
			cout << get(i,j) << " ";
		}
		cout << endl;
	}
}

template <typename type> void SpacialArray2D<type>::fill(const type& value)
{
	for (int i=0; i<=size.get_dir1()+1; ++i) {
	for (int j=0; j<=size.get_dir2()+1; ++j)
	{
		//~ cout << "in the fill loop" << endl;
		write(i,j, value);
	}}
}


template <typename type> StructuredLocalField2D<type>
	SpacialArray2D<type>::getLocalField(int i, int j) const
{
	if (i < 1 || i > size.get_dir1() || j < 1 || j > size.get_dir2())
	{
		cout << "Array Index Out of range!" << endl;
		exit(1);

		return StructuredLocalField2D<type>();
	}
	else
	{
		StructuredLocalField2D<type> field;
		field.P=get(i,j);
		field.N=get(i,j+1);
		field.S=get(i,j-1);
		field.E=get(i+1,j);
		field.W=get(i-1,j);

		return field;
	}
}

template <typename type> void SpacialArray2D<type>::set_NeumannBdyValues_all()
{
	for (int i=1; i<=size.get_dir1(); ++i)
	{
		int j=0;
		write( i,j, get(i,j+1) );

		j=size.get_dir2()+1;
		write( i,j, get(i,j-1) );
	}

	for (int j=1; j<=size.get_dir2(); ++j)
	{
		int i=0;
		write( i,j, get(i+1,j) );

		i=size.get_dir1()+1;
		write( i,j, get(i-1,j) );
	}
	
}

template <typename type> void SpacialArray2D<type>::set_DirichletBdyValues_all(const type& value)
{
	for (int i=1; i<=size.get_dir1(); ++i)
	{
		int j=0;
		write( i,j, value );

		j=size.get_dir2()+1;
		write( i,j, value );
	}

	for (int j=1; j<=size.get_dir2(); ++j)
	{
		int i=0;
		write( i,j, value );

		i=size.get_dir1()+1;
		write( i,j, value );
	}
}

template <typename type> type SpacialArray2D<type>::get_max() const	
{
	type maxValue=get(1,1);
	for (int i=1; i<=size.get_dir1(); ++i)
	for (int j=1; j<=size.get_dir2(); ++j)
		{
			type localValue=get(i,j);
			if (localValue > maxValue)
				maxValue = localValue;
		}
	return maxValue;
}

template <typename type> void SpacialArray2D<type>::vtkOutput(ofstream& output, const char* name) const
{
	output << "\nSCALARS " << name << " float" << endl;
	output << "LOOKUP_TABLE default" << endl;
	for (int j=1; j<=getSize_dir2(); ++j)
	for (int i=1; i<=getSize_dir1(); ++i)
	{
		output << get(i,j) << endl;
	}
}

template<> void SpacialArray2D<Vector2D<double> >::vtkOutput(ofstream& output, const char* name) const
{
	output << "\nVECTORS " << name << " float" << endl;
	for (int j=1; j<=getSize_dir2(); ++j)
	for (int i=1; i<=getSize_dir1(); ++i)
	{
		output 	<< get(i,j).get_dir1() << "\t"
				<< get(i,j).get_dir2() << "\t 0.0" << endl;
	}
}


template class SpacialArray2D<Vector2D<double> >;
template class SpacialArray2D<double>;
template class SpacialArray2D<int>;



