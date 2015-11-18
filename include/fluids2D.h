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
		double gamma=1.4,
		double particleMass=28.97/(1000.0*c_Avagadro),
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
			kinematicViscosity(size),
			gamma(gamma),
			particleMass(particleMass),
			particleRadius(particleRadius) {}



	//***********************************************************************
	//Print Data*************************************************************
	inline void print_molarDensity() const { molarDensity.print(); }
	inline void print_massDensity() const { massDensity.print(); }
	inline void print_temperature() const { temperature.print(); }
	inline void print_pressure() const { pressure.print(); }
	inline void print_thermalConductivity() const { thermalConductivity.print(); }

	//***********************************************************************
	//Fill Data**************************************************************
	inline void fill_temperature(const double& scalar) { temperature.fill(scalar); }
	inline void fill_velocity(const Vector2D<double>& vector) { velocity.fill(vector); }
	inline void fill_pressure(const double& scalar) { pressure.fill(scalar); }

	void init_fromBasicProps()
	{
		for (int i=0; i<=size.get_dir0()+1; ++i)
		for (int j=0; j<=size.get_dir1()+1; ++j)
		{
			double localTemp=temperature.get(i,j); 	//units: K
			double localPress=pressure.get(i,j);	//units: Pa
			Vector2D<double> localVelocity=velocity.get(i,j);	//units: m/s


			double localMolarDensity=localPress/(c_univGasConst*localTemp);
				//units: mol/m^3

			double localMassDensity=localMolarDensity*c_Avagadro*particleMass;
				//units: kg/m^3
				
			Vector2D<double> localMomentum=localMassDensity*localVelocity;
				//units: [kg/(m^2 s)]
			
			double localInternalEnergy=localMolarDensity*Cv_J_per_mol_K()
				*localTemp; //units: J/m^3

			double localKineticEnergy=0.5*localMassDensity
				*pow(localVelocity.get_magnitude(), 2.0); //units: J/m^3

			//Dependent Variables********************************************
			molarDensity.write(i,j, localMolarDensity);
			massDensity.write(i,j, localMassDensity);
			momentum.write(i,j, localMomentum);
			internalEnergy.write(i,j, localInternalEnergy);
			totalEnergy.write(i,j, localInternalEnergy+localKineticEnergy);


			//Right Hand Sides***********************************************
			continuity_rhs.write(i,j, 0.0);
			momentum_rhs.write(i,j, Vector2D<double>(0.0,0.0));
			totalEnergy_rhs.write(i,j, 0.0);

			//Properties*****************************************************
			soundSpeed.write(i,j, getSoundSpeed(i,j));
			thermalConductivity.write(i,j, getThermalConductivity(i,j));
			thermalDiffusivity.write(i,j, getThermalDiffusivity(i,j));
			dynamicViscosity.write(i,j, getDynamicViscosity(i,j));
			kinematicViscosity.write(i,j, getKinematicViscosity(i,j));
		}
	}

	//***********************************************************************
	//Get Data***************************************************************
	inline double getPressure(int i, int j) const //units: Pa
	{
		return molarDensity.get(i,j)*c_univGasConst*temperature.get(i,j);
			//units: Pa
	}
	inline double getInternalEnergy(int i, int j) const //units: J/m^3
	{
		return molarDensity.get(i,j)*Cv_J_per_mol_K()*temperature.get(i,j);
			//units: J/m^3
	}
	inline Vector2D<int> getSize() const { return size; }

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
	inline double collFreq(int i, int j) const
	{
		//Woods (1993) (pg. 40)
		return 1.414213562*c_pi*molarDensity.get(i,j)*c_Avagadro
			*thermalVelocity_avg(i,j)
			*4.0*particleRadius*particleRadius; // units: 1/s
	}


	//***********************************************************************
	//Properties*************************************************************
	/*
	 * 	Function: getThermalConductivity()
	 *	Purpose: Return the thermal conductivity based on
	 * 		self collisions only and assuming electrically neutral.
	 */
	inline double getThermalConductivity(int i, int j) const
	{
		//See Bukowski (1996).
		//See Woods (1993) pg. 64.
		return (5.0/2.0)*molarDensity.get(i,j)*c_Avagadro
			*c_k*c_k*temperature.get(i,j)
			/(particleMass*collFreq(i,j)); //units: [W m^{-1} K{-1}]
	}

	/*
	 * 	Function: getThermalDiffusivity()
	 *	Purpose: Return the thermal diffusivity based on
	 * 		the local thermal conductivity and density.
	 */
	inline double getThermalDiffusivity(int i, int j) const //units: m^2/s
	{
		return getThermalConductivity(i,j)
			/(massDensity.get(i,j)*Cv_J_per_kg_K());
			//units: m^2/s
	}

	/*
	 * 	Function: getDynamicViscosity()
	 *	Purpose: Return the dynamic viscosity based on
	 * 		the assumption of electrical neutrality.
	 */
	inline double getDynamicViscosity(int i, int j) const //units: [Pa s]
	{
		//See Woods (1993) pg. 50
		return 1.4962*(2.0/3.0)
			*(1.0/(c_pi*4.0*pow(particleRadius,2.0)))
			*sqrt( particleMass*c_k*temperature.get(i,j)/c_pi);
			//units: Pa s
	}

	/*
	 * 	Function: getKinematicViscosity()
	 *	Purpose: Return the kinematic viscosity based on
	 * 		the dynamic viscosity and the density.
	 */
	inline double getKinematicViscosity(int i, int j) const //units: [m^2/s]
	{
		//See Woods (1993) pg. 50
		return getDynamicViscosity(i,j)
			/(massDensity.get(i,j));
			//units: [m^2/s]
	}


	/*
	 * Function: getSoundSpeed()
	 * Purpose: Return the local fluid sound speed.
	 */
	inline double getSoundSpeed(int i, int j) const
	{
		return sqrt(gamma*c_k*temperature.get(i,j)/particleMass);
			// units: m/s
	}
	
	/*
	 * Function: Cv_J_per_kg_K()
	 * Purpose: Return the specific heat of the fluid at a constant volume.
	 * 	Units are J/(kg K).
	 */
	inline double Cv_J_per_kg_K() const //units: J/(kg K)
	{
		return c_k/((gamma-1.0)*particleMass); //units: J/(kg K)
	}

	/*
	 * Function: Cv_J_per_mol_K()
	 * Purpose: Return the specific heat of the fluid at a constant volume.
	 * 	Units are J/(mol K).
	 */
	inline double Cv_J_per_mol_K() const //units: J/(mol K)
	{
		return c_univGasConst/(gamma-1.0); //units: J/(mol K)
	}

};


class ThreeComponentPlasma2D
{
	private:
	Vector2D<int> size;
	
	PerfectGas2D elec;
	PerfectGas2D ions;
	PerfectGas2D neut;

	public:
	ThreeComponentPlasma2D(
		Vector2D<int> size=Vector2D<int>(10,10),
		double gamma=5.0/3.0,
		double particleMass=2.3258671e-26,
		double particleRadius=65.0e-12)
		:	elec(size, 5.0/3.0, c_eMass, 0.0),
			ions(size, gamma, particleMass, particleRadius),
			neut(size, gamma, particleMass, particleRadius) {}
};


#endif // INCLUDE_FLUIDS2D_H_
