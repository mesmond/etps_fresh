//***************************************************************************
//***************************************************************************
/*
 * File: fluids2D.h
 * Description: This file contains declarations and definitions of
 * 	several useful functions and classes related to fluids.  This includes
 * 	relevant spacial arrays and operations.
*/
//***************************************************************************
//***************************************************************************
#ifndef INCLUDE_FLUIDS2D_H_
#define INCLUDE_FLUIDS2D_H_

#include <iostream>

#include "dataTypes2D.h"
#include "simulationConstants.h"

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
	SpacialArray2D<double> soundSpeed;			//m/s
	SpacialArray2D<double> thermalConductivity;	//W/(m K)
	SpacialArray2D<double> thermalDiffusivity;	//m^2/s
	SpacialArray2D<double> dynamicViscosity;	//Pa s
	SpacialArray2D<double> kinematicViscosity;	//m^2/s


	double gamma;						//Specific Heat Ratio (constant)
	double particleMass;				//kg/particle
	double particleRadius;				//m
	double molarWeight;					//kg/mol

	double maxSoundSpeed;				//m/s

	public:
	PerfectGas2D(
		Vector2D<int> size=Vector2D<int>(10,10),
		double gamma=5.0/3.0,
		double particleMass=2.3258671e-26,
		double particleRadius=65.0e-12)
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
			thermalConductivity(size),
			thermalDiffusivity(size),
			dynamicViscosity(size),
			gamma(gamma),
			particleMass(particleMass),
			particleRadius(particleRadius) {}

	inline Vector2D<int> getSize() const { return size; }

	inline void print_molarDensity() const { molarDensity.print(); }
	inline void print_massDensity() const { massDensity.print(); }

	//***********************************************************************
	//Thermal Speeds*********************************************************
	/*
	 * 	Function: thermalVelocity_avg()
	 * 	Purpose: Return the maxwellian averaged thermal velocity.
	 */
	inline double thermalVelocity_avg(int i, int j) const
		{ return sqrt(8.0*c_k*temperature.get(i,j)/(c_pi*particleMass) ); }

	/*
	 * 	Function: thermalVelocity()
	 * 	Purpose: Return the maxwellian thermal velocity.  This is the
	 * 		most probable speed of an individual particle in the gas
	 */
	inline double thermalVelocity(int i, int j) const
		{ return sqrt(2.0*c_k*temperature.get(i,j)/particleMass ); }

	/*
	 * 	Function: thermalVelocity_rms()
	 * 	Purpose: Return the RMS value of the velocity for a
	 * 		maxwellian distribution.
	 */
	inline double thermalVelocity_rms(int i, int j) const
		{ return sqrt(3.0*c_k*temperature.get(i,j)/particleMass ); }


	//***********************************************************************
	//Collisions*************************************************************

	/*
	 * 	Function: collFreq()
	 * 	Purpose: Return the self collision frequency for the fluid assuming
	 * 		electrically neutral.
	 */
	inline double collFreq(int i, int j)
	{
		//Woods (1993) (pg. 40)
		return 1.414213562*c_pi*molarDensity.get(i,j)*c_Avagadro
			*thermalVelocity_avg(i,j)
			*4.0*particleRadius*particleRadius;
	}


	//***********************************************************************
	//Derived Properties*****************************************************
	/*
	 * 	Function: getThermalConductivty()
	 *	Purpose: Return the thermal conductivity based on
	 * 		self collisions only and assuming electrically neutral.
	 */
	inline double getThermalConductivty(int i, int j)
	{
		//See Bukowski (1996).
		return (5.0/2.0)*molarDensity.get(i,j)*c_Avagadro
			*c_k*c_k*temperature.get(i,j)
			/(particleMass*collFreq(i,j)); // [W m^{-1} K{-1}]
	}

	/*
	 * Function: getSoundSpeed()
	 * Purpose: Return the local fluid sound speed.
	 */
	inline double getSoundSpeed(int i, int j)
	{
		return sqrt(gamma*c_k*temperature.get(i,j)/particleMass);
	}
	

};


//~ class ThreeComponentPlasma2D
//~ {
	//~ private:
	//~ Vector2D<int> size;
	//~ 
	//~ PerfectGas2D elec;
	//~ PerfectGas2D ions;
	//~ PerfectGas2D neut;
//~ 
	//~ public:
	//~ ThreeComponentPlasma2D(
		//~ Vector2D<int> size=Vector2D<int>(10,10),
		//~ 
		//~ double gamma=5.0/3.0,
		//~ double particleMass=2.3258671e-26,
		//~ double particleRadius=65.0e-12)
//~ }


#endif // INCLUDE_FLUIDS2D_H_
