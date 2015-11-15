//***************************************************************************
//***************************************************************************
/*
 * File: simulationConstants.h
 * Description: This file contains several useful constants.
*/
//***************************************************************************
//***************************************************************************
#ifndef INCLUDE_SIMULATIONCONSTANTS_H_
#define INCLUDE_SIMULATIONCONSTANTS_H_

//***************************************************************************
//Useful Conversions*********************************************************
//	1 kJ/mol=96.4869 eV/atom
// 	1 J = 6.24150934e18 eV
//	1 ev = 1.60217657e-19 J
// 	1 Ohm = 1 kg m^2 s^{-3} A^{-2}  =W A^{-2}
//	1eV = 11604.50520 K
//	1 Tesla = 10^4 Gauss
//	1 amu = 1.660538921e−27 kg

//***************************************************************************
//Unit Conversions***********************************************************
// 1 Newton = 1 kg m s^-2
// 1 Joule = 1 W s
// 1 Amp = 1 C s^-1
// 1 Volt = 1 kg m^2 s^-2 C^-1  = J/C
// 1 tesla = 1 kg s^-1 C^-1 
// 1 farad m^-1 = 1 s^2 C^2 kg^-1 m^-3
// 1 henry m^-1 = 1 kg m C^-2

//***************************************************************************
//Physical Constants*********************************************************
const double c_mu0=1.25663706e-6;	
	//!< Permeability of Free Space [m kg s^{-2} A^{-2}].
const double c_eps0=8.85418782e-12;	
	//!< Permittivity of Free Space [s^4 A^2 m^{-3} kg^{-1}].
const double c_eMass=9.10938291e-31;
	//!< Mass of an electron [kg].
const double c_eCh=1.60217657e-19;	
	//!< Charge of an electron [C].
const double c_pi=3.1415926535897932384626433832795L;
	//!< Ratio of a circle's circumference to its diameter.
const double c_k=1.3806488e-23;	
	//!< Boltzmann Constant [m^2 kg s^{-2} K^{-1}] or [J/K].
const double c_sigma=5.670373e-8; 	
	//!< Stefan-Boltzmann Constant [W m^{-2} K^{-4}]
const double c_Planck=6.62606957e-34;	
	//!< Planck's Constant [m^2 kg s^{-1}] or [J s].
const double c_Avagadro=6.0221413e23; 	
	//!< Avagadro's Number [particles per mole].
const double c_LightSpeed=2.99792458e8; 
	//!< Speed of light [m/s].
const double c_protonMass=1.67262178e-27; 
	//!< Proton mass [kg].
const double c_univGasConst=c_k*c_Avagadro; 
	//!< Universal Gas Constant [J mol^{-1} K^{-1}].
const double c_eMolMass=c_eMass*c_Avagadro; 
	//!< Electron Molar Weight [kg/mol].
const double c_CoulombConstant=8.9875517873681764e9;
	//!< \brief Coulomb Constant [N m^2 C^{-2}]. 
	//!< (Used for some conversion between CGS and MKS Units.)

//***************************************************************************
//Conversion Factors*********************************************************
const double c_kg_per_amu=1.660538921e-27;
	//!< Conversion Factor [kg/amu].

//***************************************************************************
//Cell Normal Vectors********************************************************
const double cellNormal_N=1.0L;
	//!< Cell north face normal vector.
const double cellNormal_S=-1.0L;
	//!< Cell south face normal vector.
const double cellNormal_E=1.0L;
	//!< Cell east face normal vector.
const double cellNormal_W=-1.0L;
	//!< Cell west face normal vector.

#endif  // INCLUDE_SIMULATIONCONSTANTS_H_


