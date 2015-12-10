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
#include <string.h>

#include "simulationConstants.h"
#include "dataTypes2D.h"
#include "structuredGeometry2D.h"
#include "mesmond-utils.h"

class Euler2D;

using namespace std;

class PerfectGas2D
{
	friend class Euler2D;

	private:
	StructuredGeometry2D* geometry;

	//Dependent Variables****************************************************
	SpacialArray2D<double> molarDensity; 					//mol/m^3
	SpacialArray2D<double> massDensity;						//kg/m^3
	
	SpacialArray2D<Vector2D<double> > velocity;				//m/s
	SpacialArray2D<Vector2D<double> > momentum;				//kg/(m^2 s)
	
	SpacialArray2D<double> temperature;						//Kelvin
	SpacialArray2D<double> internalEnergy;					//J/m^3
	SpacialArray2D<double> totalEnergy;						//J/m^3

	SpacialArray2D<double> pressure;						//Pa

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
	//***********************************************************************
	//Constructors***********************************************************
	PerfectGas2D(
		StructuredGeometry2D& geom,
		double gamma=1.4,
		double particleMass=28.97/(1000.0*c_Avagadro),
		double particleRadius=65.0e-12); //Constructor
	PerfectGas2D(const PerfectGas2D& that); //Copy Constructor
	PerfectGas2D operator=(const PerfectGas2D& rhs) = delete; //Copy Assign
	PerfectGas2D(PerfectGas2D&& other) = delete; //Move Con
	PerfectGas2D& operator=(PerfectGas2D&& other) = delete; //Move Assign
	~PerfectGas2D(); //Deconstructor

	//***********************************************************************
	//Print Data*************************************************************
	inline void print_molarDensity() const { molarDensity.print(); }
	inline void print_massDensity() const { massDensity.print(); }
	inline void print_temperature() const { temperature.print(); }
	inline void print_pressure() const { pressure.print(); }
	inline void print_thermalConductivity() const { thermalConductivity.print(); }
	inline void print_geometry() const { geometry->print(); }

	//***********************************************************************
	//Output Data************************************************************
	void vtkOutput(const char* prefix, int& outputCount) const;
	void sodProblemOutput(double simTime);

	//***********************************************************************
	//Fill Data**************************************************************
	inline void fill_temperature(const double& scalar) { temperature.fill(scalar); }
	inline void fill_pressure(const double& scalar) { pressure.fill(scalar); }
	inline void fill_velocity(const Vector2D<double>& vector) { velocity.fill(vector); }

	void init_from_temperature_pressure_velocity();
	void init_test_sedov_lowIntensity();
	void init_test_sod_shockTube();

	//***********************************************************************
	//Update Data************************************************************
	inline void continuity_write_rhs(int i, int j, double value)
		{ continuity_rhs.write(i,j, value); }
	inline void momentum_write_rhs(int i, int j, Vector2D<double> value)
		{ momentum_rhs.write(i,j, value); }
	inline void totalEnergy_write_rhs(int i, int j, double value)
		{ totalEnergy_rhs.write(i,j, value); }

	void update_from_rhs(int i, int j, double timeStep);

	void update_boundary_values();

	//***********************************************************************
	//Get Data***************************************************************
	inline StructuredLocalField2D<double> getCellAreas(int i, int j) const
		{ return geometry->getCellAreas(i,j); }

	inline double getVolume(int i, int j) const
		{ return geometry->getVolume(i,j); }

	inline Vector2D<int> getSize() const { return geometry->getSize(); }

	double get_explicit_timeStep(double CFL_number=0.5) const;


	//***********************************************************************
	//Calculate Data*********************************************************
	/*
	 * Functions:
	 * 		calcPressure() 				//units: Pa
	 * 		calcSoundSpeed()			//units: m/s
	 * 		calcThermalConductivity()	//units: [W m^{-1} K{-1}]
	 *		calcThermalDiffusivity()	//units: m^2/s
	 *		calcDynamicViscosity()		//units: [Pa s]
	 *		calcKinematicViscosity()	//units: m^2/s
	 * 
	 * Purpose: Return specified fluid properties based on
	 * 	elementary fluid values and dependent variables.
	 *
	 * Descriptions:
	 * 		calcPressure(): Returns the pressure based on the density
	 * 			and temperature.
	 * 		calcSoundSpeed(): Returns the sound speed based on the fluid
	 * 			temperature and fluid particle mass.
	 * 		calcThermalConductivity(): Returns the thermal conductivity
	 * 			based on the work of Bukowski (1996) and Woods (1993, pg. 64).
	 * 			This is based on the fluid temperature, density, and particle
	 * 			mass.
	 *		calcThermalDiffusivity(): Returns the thermal diffusivity
	 * 			based on the calculated thermal conductivity, fluid density,
	 * 			and specific heat.
	 *		calcDynamicViscosity(): Returns an approximation of the fluid's
	 * 			dynamic viscosity based on the work of Woods (1993, pg. 50).
	 *		calcKinematicViscosity(): Returns the kinematic viscosity
	 * 			based on a calculation of the dynamic viscosity and the
	 * 			the fluid density.
	 */

	double calcPressure(int i, int j) const; //units: Pa
	inline double calcSoundSpeed(int i, int j) const //units: m/s
		{ return sqrt(gamma*c_k*temperature.get(i,j)/particleMass); }
	double calcThermalConductivity(int i, int j) const; //units: [W m^{-1} K{-1}]
	double calcThermalDiffusivity(int i, int j) const; //units: m^2/s
	double calcDynamicViscosity(int i, int j) const; //units: [Pa s]
	double calcKinematicViscosity(int i, int j) const; //units: [m^2/s]


	//***********************************************************************
	//Thermal Speeds*********************************************************
	/*
	 * Functions:
	 * 		thermalVelocity_avg()		//units: m/s
	 * 		thermalVelocity()			//units: m/s
	 * 		thermalVelocity_rms()		//units: m/s
	 *
	 * Purpose: Return different representations of the velocity of
	 * 	particles with a Maxwellian distribution.
	 *
	 * Descriptions:
	 * 		thermalVelocity_avg(): Returns the average speed of all
	 * 			the fluid particles in a Maxwellian distibution.
	 * 		thermalVelocity(): Returns the most probable speed of
	 * 			a fluid particle with a Maxwellian distribution.
	 * 		thermalVelocity_rms(): Returns the RMS value of the velocity of
	 * 			fluid particles in a Maxwellian Distribution.
	 */

	inline double thermalVelocity_avg(int i, int j) const //units: m/s
		{ return sqrt(8.0*c_k*temperature.get(i,j)/(c_pi*particleMass) ); }
		
	inline double thermalVelocity(int i, int j) const //units: m/s
		{ return sqrt(2.0*c_k*temperature.get(i,j)/particleMass ); }
		
	inline double thermalVelocity_rms(int i, int j) const //units: m/s
		{ return sqrt(3.0*c_k*temperature.get(i,j)/particleMass ); }


	//***********************************************************************
	//Collisions*************************************************************
	/*
	 * 	Function:
	 * 		collFreq()			//units: 1/s
	 * 
	 * 	Purpose: Calculate and return values related to the collisionality of
	 * 		the fluid.
	 *
	 * 	Description:
	 * 		collFreq(): Returns the collision frequency
	 * 			due to intra-species collisions. For more details,
	 * 			See Woods (1993) pg. 40.
	 */

	inline double collFreq(int i, int j) const; //units: 1/s


	//***********************************************************************
	//Specific Heat**********************************************************
	/*
	 * Functions:
	 * 		Cv_J_per_kg_K()			//units: J/(kg K)
	 * 		Cv_J_per_mol_K()		//units: J/(mol K)
	 *
	 * Purpose: Return the fluid specific heat in desired units.
	 *
	 * Descriptions:
	 * 		Cv_J_per_kg_K(): Return the specific heat of
	 * 			the fluid at a constant volume. Units are J/(kg K).
	 * 		Cv_J_per_mol_K(): Return the specific heat of
	 * 			the fluid at a constant volume. Units are J/(mol K).
	 */

	inline double Cv_J_per_kg_K() const //units: J/(kg K)
		{ return c_k/((gamma-1.0)*particleMass); }

	inline double Cv_J_per_mol_K() const //units: J/(mol K)
		{ return c_univGasConst/(gamma-1.0); }

};

class Euler2D
{
	public:
	double get_continuity_rhs(int i, int j, const PerfectGas2D& fluid) const
	{
		double numericalFlux=get_numericalDissipationFlux(i,j, fluid, fluid.massDensity);

		return 	(-1.0)*fluid.geometry->divergence(i,j, fluid.massDensity,
					fluid.velocity)
				+(-1.0)*numericalFlux;
	}


	Vector2D<double> get_momentum_rhs(int i, int j, const PerfectGas2D& fluid) const
	{
		Vector2D<double> numericalFlux=get_numericalDissipationFlux(i,j, fluid, fluid.momentum);

		return 	(-1.0)*fluid.geometry->divergence(i,j, fluid.massDensity,
					fluid.velocity, fluid.velocity)
				-fluid.geometry->gradient(i,j, fluid.pressure)
				+(-1.0)*numericalFlux;
	}


	double get_totEnergy_rhs(int i, int j, const PerfectGas2D& fluid) const
	{
		double numericalFlux=get_numericalDissipationFlux(i,j, fluid, fluid.totalEnergy);

		return 	(-1.0)*fluid.geometry->divergence(i,j, fluid.totalEnergy,
					fluid.velocity)
				-fluid.geometry->divergence(i,j, fluid.pressure,
					fluid.velocity)
				+(-1.0)*numericalFlux;
	}

	private:
	template <typename T> T get_numericalDissipationFlux(int i, int j, const PerfectGas2D& fluid,
		SpacialArray2D<T> array) const
	{
		StructuredLocalField2D<T> localField;
		StructuredLocalField2D<double> soundSpeed;
		StructuredLocalField2D<double> maxSoundSpeed;
		StructuredLocalField2D<T> localFlux;

		localField=array.getLocalField(i,j);
		soundSpeed=fluid.soundSpeed.getLocalField(i,j);
		
		maxSoundSpeed.N=max_2arg(soundSpeed.P, soundSpeed.N);
		maxSoundSpeed.S=max_2arg(soundSpeed.P, soundSpeed.S);
		maxSoundSpeed.E=max_2arg(soundSpeed.P, soundSpeed.E);
		maxSoundSpeed.W=max_2arg(soundSpeed.P, soundSpeed.W);

		localFlux.N=-0.5*maxSoundSpeed.N*(localField.N-localField.P);
		localFlux.S=-0.5*maxSoundSpeed.S*(localField.P-localField.S);
		localFlux.E=-0.5*maxSoundSpeed.E*(localField.E-localField.P);
		localFlux.W=-0.5*maxSoundSpeed.W*(localField.P-localField.W);

		StructuredLocalField2D<double> area=fluid.geometry->getCellAreas(i,j);

		//Compute Divergence.
		T divergence=(1.0/fluid.geometry->getVolume(i,j))*(
			 cellNormal_N*localFlux.N*area.N
			+cellNormal_S*localFlux.S*area.S
			+cellNormal_E*localFlux.E*area.E
			+cellNormal_W*localFlux.W*area.W );

		return divergence;

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
