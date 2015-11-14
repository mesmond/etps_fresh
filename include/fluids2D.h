//***************************************************************************
//***************************************************************************
/*
 * File: fluids2D.h
 * Description: This file contains declarations and definitions of
 * 	several useful functions related to fluids.  This includes
 * 	relevant spacial arrays and operations.
*/
//***************************************************************************
//***************************************************************************
#ifndef INCLUDE_FLUIDS2D_H_
#define INCLUDE_FLUIDS2D_H_

#include <iostream>

#include "dataTypes2D.h"

using namespace std;

class Fluid2D
{
	private:
	Vector2D<int> count;

	public:
	SpacialArray2D<double> pressure;
	SpacialArray2D<Vector2D<double> > velocity;
	SpacialArray2D<double> temperature;

	SpacialArray2D<Vector2D<double> > momentum;
	SpacialArray2D<double> internalEnergy;
	SpacialArray2D<double> totalEnergy;

	Fluid2D(int size0=10, int size1=10) :
		pressure(size0, size1),
		velocity(size0, size1),
		temperature(size0, size1),
		momentum(size0, size1),
		internalEnergy(size0, size1),
		totalEnergy(size0, size1)
	{
		count.write_dir0(size0);
		count.write_dir1(size1);
	}

	Vector2D<int> getCount() const { return count; }
};

class IdealGas : public Fluid2D
{
	SpacialArray2D<double> massDensity;
	SpacialArray2D<double> molarDensity;

	//~ void updatePressure()
	//~ {
		//~ int count_dir0=this->getCount().get_dir0();
		//~ int count_dir1=this->getCount().get_dir1();
//~ 
		//~ 
	//~ }

	
	
};


#endif // INCLUDE_FLUIDS2D_H_
