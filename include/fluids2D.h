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
#include "operators2D.h"
#include "simulationConstants.h"

using namespace std;

class PerfectGas2D
{
	private:
	//Mesh Size**************************************************************
	StructuredGeometry2D geometry;

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
		StructuredGeometry2D geometry,
		double gamma=1.4,
		double particleMass=28.97/(1000.0*c_Avagadro),
		double particleRadius=65.0e-12)
		:
			geometry(geometry),
			//Dependent Variables********************************************
			molarDensity(geometry.getSize()),
			massDensity(geometry.getSize()),
			velocity(geometry.getSize()),
			momentum(geometry.getSize()),
			temperature(geometry.getSize()),
			internalEnergy(geometry.getSize()),
			totalEnergy(geometry.getSize()),
			pressure(geometry.getSize()),
			//Right hand sides***********************************************
			continuity_rhs(geometry.getSize()),
			momentum_rhs(geometry.getSize()),
			totalEnergy_rhs(geometry.getSize()),
			//Properties*****************************************************
			soundSpeed(geometry.getSize()),
			thermalConductivity(geometry.getSize()),
			thermalDiffusivity(geometry.getSize()),
			dynamicViscosity(geometry.getSize()),
			kinematicViscosity(geometry.getSize()),
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
	inline void print_geometry() const { geometry.print(); }

	//***********************************************************************
	//Fill Data**************************************************************
	inline void fill_temperature(const double& scalar) { temperature.fill(scalar); }
	inline void fill_velocity(const Vector2D<double>& vector) { velocity.fill(vector); }
	inline void fill_pressure(const double& scalar) { pressure.fill(scalar); }

	void init_fromBasicProps();

	//***********************************************************************
	//Calculate Data*********************************************************
	/*
	 * Functions: calcPressure(), calcInternalEnergy()
	 * Purpose: Calculate and return the specified value based on
	 * 	other properties in the fluid.
	 *
	 * Description:
	 * 	calcPressure() returns the pressure
	 * 		as computed from the molar density and the temperature.
	 * 	calcInternalEnergy() returns the internal fluid energy based on the
	 * 		molar density and the temperature.
	 *
	 * 	Units are as shown by each function.
	 */
	double calcPressure(int i, int j) const; //units: Pa
	double calcInternalEnergy(int i, int j) const; //units: J/m^3

	//***********************************************************************
	//Get Data***************************************************************
	inline Vector2D<int> getSize() const { return geometry.getSize(); }

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
	 *
	 * 	Description: See Woods (1993) pg. 40.
	 */
	inline double collFreq(int i, int j) const; //units: 1/s


	//***********************************************************************
	//Properties*************************************************************
	/*
	 * 	Function: getThermalConductivity()
	 *	Purpose: Return the thermal conductivity based on
	 * 		self collisions only and assuming electrically neutral.
	 */
	double getThermalConductivity(int i, int j) const; //units: [W m^{-1} K{-1}]

	/*
	 * 	Function: getThermalDiffusivity()
	 *	Purpose: Return the thermal diffusivity based on
	 * 		the local thermal conductivity and density.
	 */
	double getThermalDiffusivity(int i, int j) const; //units: m^2/s

	/*
	 * 	Function: getDynamicViscosity()
	 *	Purpose: Return the dynamic viscosity based on
	 * 		the assumption of electrical neutrality.
	 *
	 * 	Description: See Woods (1993) pg. 50
	 */
	double getDynamicViscosity(int i, int j) const; //units: [Pa s]

	/*
	 * 	Function: getKinematicViscosity()
	 *	Purpose: Return the kinematic viscosity based on
	 * 		the dynamic viscosity and the density.
	 */
	double getKinematicViscosity(int i, int j) const; //units: [m^2/s]


	/*
	 * Function: getSoundSpeed()
	 * Purpose: Return the local fluid sound speed.
	 */
	inline double getSoundSpeed(int i, int j) const //units: m/s
		{ return sqrt(gamma*c_k*temperature.get(i,j)/particleMass); }
	
	/*
	 * Function: Cv_J_per_kg_K()
	 * Purpose: Return the specific heat of the fluid at a constant volume.
	 * 	Units are J/(kg K).
	 */
	inline double Cv_J_per_kg_K() const //units: J/(kg K)
		{ return c_k/((gamma-1.0)*particleMass); }

	/*
	 * Function: Cv_J_per_mol_K()
	 * Purpose: Return the specific heat of the fluid at a constant volume.
	 * 	Units are J/(mol K).
	 */
	inline double Cv_J_per_mol_K() const //units: J/(mol K)
		{ return c_univGasConst/(gamma-1.0); }

};

class Euler2D
{
	public:
	Operators2D* operate;
	Euler2D(Operators2D* ptr)
	{
		operate=ptr;

		cout << "operate->getVolume(i,j)=" << operate->getVolume(2,3) << endl;
	}



	double continuity_rhs(const PerfectGas2D& fluid)
	{
		return 0.0;
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
		//~ double gamma=5.0/3.0,
		//~ double particleMass=2.3258671e-26,
		//~ double particleRadius=65.0e-12)
		//~ :	elec(size, 5.0/3.0, c_eMass, 0.0),
			//~ ions(size, gamma, particleMass, particleRadius),
			//~ neut(size, gamma, particleMass, particleRadius) {}
//~ };


#endif // INCLUDE_FLUIDS2D_H_
