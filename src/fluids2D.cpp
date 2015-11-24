//***************************************************************************
//***************************************************************************
/*
 * File: fluids2D.cpp
 * Description: This file contains class definitions and implementations
 * 	for types related to fluids.
*/
//***************************************************************************
//***************************************************************************
#include "fluids2D.h"

//***************************************************************************
//Fill Data******************************************************************
void PerfectGas2D::init_fromBasicProps()
{
	for (int i=0; i<=getSize().get_dir1()+1; ++i)
	for (int j=0; j<=getSize().get_dir2()+1; ++j)
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

//***************************************************************************
//Calculate Data*************************************************************
double PerfectGas2D::calcPressure(int i, int j) const //units: Pa
{
	return molarDensity.get(i,j)*c_univGasConst*temperature.get(i,j);
		//units: Pa
}
double PerfectGas2D::calcInternalEnergy(int i, int j) const //units: J/m^3
{
	return molarDensity.get(i,j)*Cv_J_per_mol_K()*temperature.get(i,j);
		//units: J/m^3
}

//***************************************************************************
//Collisions*****************************************************************
double PerfectGas2D::collFreq(int i, int j) const
{
	//Woods (1993) (pg. 40)
	return 1.414213562*c_pi*molarDensity.get(i,j)*c_Avagadro
		*thermalVelocity_avg(i,j)
		*4.0*particleRadius*particleRadius; // units: 1/s
}

//***************************************************************************
//Properties*****************************************************************
double PerfectGas2D::getThermalConductivity(int i, int j) const //units: [W m^{-1} K{-1}]
{
	//See Bukowski (1996).
	//See Woods (1993) pg. 64.
	return (5.0/2.0)*molarDensity.get(i,j)*c_Avagadro
		*c_k*c_k*temperature.get(i,j)
		/(particleMass*collFreq(i,j)); //units: [W m^{-1} K{-1}]
}

double PerfectGas2D::getThermalDiffusivity(int i, int j) const //units: m^2/s
{
	return getThermalConductivity(i,j)
		/(massDensity.get(i,j)*Cv_J_per_kg_K());
		//units: m^2/s
}

double PerfectGas2D::getDynamicViscosity(int i, int j) const //units: [Pa s]
{
	//See Woods (1993) pg. 50
	return 1.4962*(2.0/3.0)
		*(1.0/(c_pi*4.0*pow(particleRadius,2.0)))
		*sqrt( particleMass*c_k*temperature.get(i,j)/c_pi);
		//units: Pa s
}

double PerfectGas2D::getKinematicViscosity(int i, int j) const //units: [m^2/s]
{
	return getDynamicViscosity(i,j)
		/(massDensity.get(i,j));
		//units: [m^2/s]
}
