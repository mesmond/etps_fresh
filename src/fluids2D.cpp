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

PerfectGas2D::PerfectGas2D( //Constructor
	StructuredGeometry2D& geom,
	double gamma,
	double particleMass,
	double particleRadius)
	:	//Dependent Variables************************************************
		molarDensity(geom.getSize()),
		massDensity(geom.getSize()),
		velocity(geom.getSize()),
		momentum(geom.getSize()),
		temperature(geom.getSize()),
		internalEnergy(geom.getSize()),
		totalEnergy(geom.getSize()),
		pressure(geom.getSize()),
		//Right hand sides***************************************************
		continuity_rhs(geom.getSize()),
		momentum_rhs(geom.getSize()),
		totalEnergy_rhs(geom.getSize()),
		//Properties*********************************************************
		soundSpeed(geom.getSize()),
		thermalConductivity(geom.getSize()),
		thermalDiffusivity(geom.getSize()),
		dynamicViscosity(geom.getSize()),
		kinematicViscosity(geom.getSize()),
		gamma(gamma),
		particleMass(particleMass),
		particleRadius(particleRadius)
{
	//Check to see if the input Structured Geometry is Cylindrical.
	StructuredCylGeometry2D* test=dynamic_cast<StructuredCylGeometry2D*>(&geom);
	if (test)
	{
		cout << "Initializing PerfectGas2D in Cylindrical Coordinates." << endl;
		geometry = new StructuredCylGeometry2D(geom);
	}
	else
	{
		cout << "Initializing PerfectGas2D in Rectangular Coordinates." << endl;
		geometry = new StructuredGeometry2D(geom);
	}
}

PerfectGas2D::PerfectGas2D(const PerfectGas2D& that) //Copy Constructor
	:	//Dependent Variables************************************************
		molarDensity(that.geometry->getSize()),
		massDensity(that.geometry->getSize()),
		velocity(that.geometry->getSize()),
		momentum(that.geometry->getSize()),
		temperature(that.geometry->getSize()),
		internalEnergy(that.geometry->getSize()),
		totalEnergy(that.geometry->getSize()),
		pressure(that.geometry->getSize()),
		//Right hand sides***************************************************
		continuity_rhs(that.geometry->getSize()),
		momentum_rhs(that.geometry->getSize()),
		totalEnergy_rhs(that.geometry->getSize()),
		//Properties*********************************************************
		soundSpeed(that.geometry->getSize()),
		thermalConductivity(that.geometry->getSize()),
		thermalDiffusivity(that.geometry->getSize()),
		dynamicViscosity(that.geometry->getSize()),
		kinematicViscosity(that.geometry->getSize()),
		gamma(that.gamma),
		particleMass(that.particleMass),
		particleRadius(that.particleRadius)
{
	//Check to see if the input Structured Geometry is Cylindrical.
	StructuredCylGeometry2D* test=dynamic_cast<StructuredCylGeometry2D*>(that.geometry);
	if (test)
	{
		cout << "Copy Constructor: PerfectGas2D in Cylindrical Coordinates." << endl;
		geometry = new StructuredCylGeometry2D(*that.geometry);
	}
	else
	{
		cout << "Copy Constructor: PerfectGas2D in Rectangular Coordinates." << endl;
		geometry = new StructuredGeometry2D(*that.geometry);
	}
}

PerfectGas2D::~PerfectGas2D()
{
	if (geometry != nullptr)
	{
		delete geometry;
		geometry=nullptr;
	}
}

//***************************************************************************
//Print Data*****************************************************************
void PerfectGas2D::vtkOutput(const char* prefix, int& outputCount) const
{	
	cout << "\t\t Writing Output File, prefix=" << prefix
		<< ", index=" << outputCount << endl;

	//Prepare the output file extension.
	const char *extension=".vtk";

	//Convert the output index to a six digit integer.
	char outputCount_char[7];
	char buffer[100];

	sprintf(outputCount_char,"%.6d",outputCount);

	//Store the file name inside the buffer.
	//Start with the prefix.
	int charLocation=0;
	while (prefix[charLocation] != '\0')
	{
		buffer[charLocation]=prefix[charLocation];
		charLocation++;
	}

	//Add the pre-index notation.
	buffer[charLocation]='@';
	charLocation++;
	buffer[charLocation]='-';
	charLocation++;
	int prefixLength=charLocation;

	//Add the output index.
	while ( (charLocation-prefixLength) < 6 )
	{
		buffer[charLocation]=outputCount_char[charLocation-prefixLength];
		charLocation++;
	}

	//Add the extension.
	prefixLength=charLocation;
	while ( (charLocation-prefixLength) < 4 )
	{
		buffer[charLocation]=extension[charLocation-prefixLength];
		charLocation++;
	}

	//Add the terminating character.
	buffer[charLocation]='\0';

	//Open the file.
	ofstream output;
	output.open(buffer, ios::out); //Clear contents if it exists.

	//Output the data to the output file.
	if (output.is_open())
	{
		output << "# vtk DataFile Version 2.0" << endl;
		output 	<< "Blank"
				<< scientific
				<< setprecision(10)
				<< endl;
		output << "ASCII" << endl;

		geometry->vtkOutput(output);

		output << "\nCELL_DATA "
			<< geometry->getCount_dir1()*geometry->getCount_dir2() << endl;
		massDensity.vtkOutput(output, "massDensity");
		molarDensity.vtkOutput(output, "molarDensity");
		velocity.vtkOutput(output, "velocity");
		temperature.vtkOutput(output, "temperature");
		pressure.vtkOutput(output, "pressure");

		momentum.vtkOutput(output, "momentum");
		totalEnergy.vtkOutput(output, "totalEnergy");

		output.close();
	}
	else
	{
		cout << "Unable to open output file=" << buffer << endl;
		exit(1);
	}

	outputCount++;
}

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
		soundSpeed.write(i,j, calcSoundSpeed(i,j));
		thermalConductivity.write(i,j, calcThermalConductivity(i,j));
		thermalDiffusivity.write(i,j, calcThermalDiffusivity(i,j));
		dynamicViscosity.write(i,j, calcDynamicViscosity(i,j));
		kinematicViscosity.write(i,j, calcKinematicViscosity(i,j));
	}
}

//***************************************************************************
//Update Data****************************************************************
void PerfectGas2D::update_from_rhs(int i, int j, double timeStep)
{
	//Get new Values*********************************************************
	double massDensity_new=massDensity.get(i,j)
		+timeStep*continuity_rhs.get(i,j);
	Vector2D<double> momentum_new=momentum.get(i,j)
		+timeStep*momentum_rhs.get(i,j);
	double totalEnergy_new=totalEnergy.get(i,j)
		+timeStep*totalEnergy_rhs.get(i,j);

	//Update Dependent Variables*********************************************
	massDensity.write(i,j, massDensity_new);
	momentum.write(i,j, momentum_new);
	totalEnergy.write(i,j, totalEnergy_new);

	//Update Object State****************************************************
	molarDensity.write(i,j, massDensity_new/(particleMass*c_Avagadro) );
	velocity.write(i,j, momentum_new*(1.0/massDensity_new));
	double localKineticEnergy=0.5*massDensity_new*pow(
		velocity.get(i,j).get_magnitude(), 2.0);
	internalEnergy.write(i,j, totalEnergy_new-localKineticEnergy);
	temperature.write(i,j, internalEnergy.get(i,j)
		*(1.0/(massDensity_new*Cv_J_per_kg_K())));

	pressure.write(i,j, calcPressure(i,j));
	soundSpeed.write(i,j, calcSoundSpeed(i,j) );
	thermalConductivity.write(i,j, calcThermalConductivity(i,j) );
	thermalDiffusivity.write(i,j, calcThermalDiffusivity(i,j) );
	dynamicViscosity.write(i,j, calcDynamicViscosity(i,j) );
	kinematicViscosity.write(i,j, calcKinematicViscosity(i,j) );
}

void PerfectGas2D::update_boundary_values()
{
	massDensity.set_NeumannBdyValues_all();
	molarDensity.set_NeumannBdyValues_all();
	momentum.set_DirichletBdyValues_all(Vector2D<double>(0.0,0.0) );
	velocity.set_DirichletBdyValues_all(Vector2D<double>(0.0,0.0) );
	totalEnergy.set_NeumannBdyValues_all();
	internalEnergy.set_NeumannBdyValues_all();

	temperature.set_NeumannBdyValues_all();
}

//***************************************************************************
//Get Data*******************************************************************

double PerfectGas2D::get_explicit_timeStep(double CFL_number) const
{
	double maxSoundSpeed=soundSpeed.get_max();
	double minMeshSpacing=geometry->get_minMeshSpacing();

	return CFL_number*minMeshSpacing/maxSoundSpeed;
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
double PerfectGas2D::calcThermalConductivity(int i, int j) const //units: [W m^{-1} K{-1}]
{
	//See Bukowski (1996).
	//See Woods (1993) pg. 64.
	return (5.0/2.0)*molarDensity.get(i,j)*c_Avagadro
		*c_k*c_k*temperature.get(i,j)
		/(particleMass*collFreq(i,j)); //units: [W m^{-1} K{-1}]
}

double PerfectGas2D::calcThermalDiffusivity(int i, int j) const //units: m^2/s
{
	return calcThermalConductivity(i,j)
		/(massDensity.get(i,j)*Cv_J_per_kg_K());
		//units: m^2/s
}

double PerfectGas2D::calcDynamicViscosity(int i, int j) const //units: [Pa s]
{
	//See Woods (1993) pg. 50
	return 1.4962*(2.0/3.0)
		*(1.0/(c_pi*4.0*pow(particleRadius,2.0)))
		*sqrt( particleMass*c_k*temperature.get(i,j)/c_pi);
		//units: Pa s
}

double PerfectGas2D::calcKinematicViscosity(int i, int j) const //units: [m^2/s]
{
	return calcDynamicViscosity(i,j)
		/(massDensity.get(i,j));
		//units: [m^2/s]
}

//***************************************************************************
//Thermal Speeds*************************************************************


//***************************************************************************
//Collisions*****************************************************************


//***************************************************************************
//Properties*****************************************************************
