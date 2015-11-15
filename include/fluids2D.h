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

class PerfectGas2D
{
	private:
	//Mesh Size**************************************************************
	Vector2D<int> size;

	//Dependent Variables****************************************************
	SpacialArray2D<double> molarDensity; 					//mol/m^3
	SpacialArray2D<double> massDensity;						//kg/m^3
	
	SpacialArray2D<Vector2D<double> > velocity;				//m/s
	SpacialArray2D<Vector2D<double> > momentum;				//kg/(m^2 s)
	
	SpacialArray2D<double> temperature;						//Kelvin
	SpacialArray2D<double> internalEnergy;					//J/m^3
	SpacialArray2D<double> totalEnergy;						//J/m^3
	
	SpacialArray2D<double> pressure;

	//Governing Equations' right hand sides**********************************
	SpacialArray2D<double> continuity_rhs;
	SpacialArray2D<Vector2D<double> > momentum_rhs;
	SpacialArray2D<double> totalEnergy_rhs;

	//Properties*************************************************************
	SpacialArray2D<double> soundSpeed;	//m/s

	double gamma;						//Specific Heat Ratio (constant)
	double particleMass;				//kg
	double maxSoundSpeed;				//m/s

	public:
	PerfectGas2D(
		Vector2D<int> size=Vector2D<int>(10,10),
		double particleMass=2.66e-26, double gamma=1.6667)
		:
			size(size),
			//Dependent Variables********************************************
			molarDensity(size),
			massDensity(size),
			velocity(size),
			momentum(size),
			temperature(size),
			internalEnergy(size),
			totalEnergy(size),
			pressure(size),
			//Right hand sides***********************************************
			continuity_rhs(size),
			momentum_rhs(size),
			totalEnergy_rhs(size),
			//Properties*****************************************************
			soundSpeed(size),
			gamma(gamma),
			particleMass(particleMass) {}

	void print_massDensity() const { massDensity.print(); }
};


#endif // INCLUDE_FLUIDS2D_H_
