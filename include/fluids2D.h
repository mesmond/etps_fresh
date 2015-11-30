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
//~ #include <typeinfo>

#include "simulationConstants.h"
#include "dataTypes2D.h"
#include "structuredGeometry2D.h"

class Euler2D;

using namespace std;

class PerfectGas2D
{
	friend class Euler2D;
	
	private:
	//Mesh Size**************************************************************
	StructuredGeometry2D* geometry;

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
	//Fill Data**************************************************************
	inline void fill_temperature(const double& scalar) { temperature.fill(scalar); }
	inline void fill_velocity(const Vector2D<double>& vector) { velocity.fill(vector); }
	inline void fill_pressure(const double& scalar) { pressure.fill(scalar); }

	void init_fromBasicProps();



	void init_test()
	{
		double air_temperature=20.0;	//deg C
		fill_temperature(air_temperature+273.15);	//K
		fill_pressure(101325.0);					//Pa
		fill_velocity(Vector2D<double>(0.0,0.0));	//m/s
		
		init_fromBasicProps();

		temperature.write(1,1, 800);
		pressure.write(1,1, calcPressure(1,1));
	}

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


	/*
	 * Function: calcSoundSpeed()
	 * Purpose: Return the local fluid sound speed based on the
	 * 	temperature, specific heat ratio, and particle mass.
	 */
	inline double calcSoundSpeed(int i, int j) const //units: m/s
		{ return sqrt(gamma*c_k*temperature.get(i,j)/particleMass); }

	/*
	 * 	Function: calcThermalConductivity()
	 *	Purpose: Return the thermal conductivity based on
	 * 		self collisions only and assuming electrically neutral.
	 */
	double calcThermalConductivity(int i, int j) const; //units: [W m^{-1} K{-1}]

	/*
	 * 	Function: calcThermalDiffusivity()
	 *	Purpose: Return the thermal diffusivity based on
	 * 		the local thermal conductivity and density.
	 */
	double calcThermalDiffusivity(int i, int j) const; //units: m^2/s

	/*
	 * 	Function: calcDynamicViscosity()
	 *	Purpose: Return the dynamic viscosity based on
	 * 		the assumption of electrical neutrality.
	 *
	 * 	Description: See Woods (1993) pg. 50
	 */
	double calcDynamicViscosity(int i, int j) const; //units: [Pa s]

	/*
	 * 	Function: calcKinematicViscosity()
	 *	Purpose: Return the kinematic viscosity based on
	 * 		the dynamic viscosity and the density.
	 */
	double calcKinematicViscosity(int i, int j) const; //units: [m^2/s]


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
