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
		particleRadius(particleRadius),
		molarWeight(particleMass*c_Avagadro),
		maxSoundSpeed(0.0)
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
		particleRadius(that.particleRadius),
		molarWeight(particleMass*c_Avagadro),
		maxSoundSpeed(0.0)
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

//***************************************************************************
//Output Data****************************************************************
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

void PerfectGas2D::sodProblemOutput(double simTime)
{
	char zoneCount_char[9];
	sprintf(zoneCount_char,"%.8d", getSize().get_dir1());
	string prefix="SodProblemOutput-";
	string extension=".out";
	string number=zoneCount_char;
	string fileName="SodProblemOutput-"+number+"z"+".out";
	
	//Open the file.
	ofstream output;
	output.open(fileName.c_str(), ios::out); //Clear contents if it exists.

	//Output the data to the output file.
	if (output.is_open())
	{
		output 	<< scientific
				<< setprecision(10);
		
		output << "# Sod Problem Data" << endl;
		output << "# Designed for a Shock propagating in the dir1 direction." << endl;
		output << "# Number of zones in direction of interest="
			<< getSize().get_dir1() << endl;
		output << "# Output Time(s)=" << simTime << endl;

		//Column Headers
		output 	<< "#Position_dir1(cm)"
				<< "\t Pressure(Mbar)"
				<< "\t MassDensity(g/cc)"
				<< "\t Veclocity_dir1(cm/micsec)"
				<< "\t InternalEnergy(Mbar cc/g)"
				<< endl;

		int j=1;
		for (int i=1; i<=getSize().get_dir1(); ++i)
		{
			output 	<< geometry->getPoint(i,j).get_dir1()*1.0e2 //units: cm
					<< "\t" << pressure.get(i,j)*1.0e-11 //units: Mbar
					<< "\t" << massDensity.get(i,j)*1.0e-3 //units: g/cc
					<< "\t" << velocity.get(i,j).get_dir1()*1.0e-4 //units: cm/micsec
					<< "\t" << temperature.get(i,j)*Cv_J_per_kg_K()*1.0e-8 //units: Mbar cc/g
					<< endl;
		}

		output.close();
	}
	else
	{
		cout << "Unable to open output file=" << fileName << endl;
		exit(1);
	}
}



//***************************************************************************
//Fill Data******************************************************************
void PerfectGas2D::init_from_temperature_pressure_velocity()
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

void PerfectGas2D::init_test_sedov_lowIntensity()
{
	double air_temperature=20.0;	//deg C
	fill_temperature(air_temperature+273.15);	//K
	fill_pressure(101325.0);					//Pa
	fill_velocity(Vector2D<double>(0.0,0.0));	//m/s

	init_from_temperature_pressure_velocity();

	temperature.write(1,1, 800);
	pressure.write(1,1, calcPressure(1,1));

	init_from_temperature_pressure_velocity();
}

void PerfectGas2D::init_test_sod_shockTube()
{
	fill_velocity(Vector2D<double>(0.0,0.0));	//m/s

	double highPressure=1.0e11; 	//units: Pa
	double highTemp=92676.71837;	//units: K

	double lowPressure=0.1e11; 		//units: Pa
	double lowTemp=74141.37469; 	//units: K

	fill_temperature(lowTemp);	//K
	fill_pressure(lowPressure);					//Pa

	init_from_temperature_pressure_velocity();

	for (int i=0; i<=getSize().get_dir1()+1; ++i)
	for (int j=0; j<=getSize().get_dir2()+1; ++j)
	{
		if (i<=getSize().get_dir1()/2)
		{
			pressure.write(i,j, highPressure);
			temperature.write(i,j, highTemp);
		}
	}

	init_from_temperature_pressure_velocity();
}




//***************************************************************************
//Update Data****************************************************************
void PerfectGas2D::update_from_rhs(int i, int j, double timeStep)
{

	//~ if (i == getSize().get_dir1()/2 && j==1)
	//~ {
		//~ cout << "MassDensity Before=" << massDensity.get(i,j);
		//~ cout << "\t MassDensity RHS=" << continuity_rhs.get(i,j);
		//~ cout << "\t (i=" << i << ")" << endl;
	//~ }


	
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


	//~ if (i == getSize().get_dir1()/2 && j==1)
	//~ {
		//~ cout << "MassDensity After=" << massDensity.get(i,j);
		//~ cout << "\t MassDensity RHS=" << continuity_rhs.get(i,j);
		//~ cout << "\t (i=" << i << ")" << endl;
//~ 
		//~ cin.get();
	//~ }



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
	//~ momentum.set_DirichletBdyValues_all(Vector2D<double>(0.0,0.0) );
	momentum.set_NeumannBdyValues_all();
	velocity.set_DirichletBdyValues_all(Vector2D<double>(0.0,0.0) );
	totalEnergy.set_NeumannBdyValues_all();
	internalEnergy.set_NeumannBdyValues_all();

	temperature.set_NeumannBdyValues_all();
	pressure.set_NeumannBdyValues_all();
}

//***************************************************************************
//Get Data*******************************************************************

double PerfectGas2D::get_explicit_timeStep(double CFL_number, double compression_number)
{
	maxSoundSpeed=soundSpeed.get_max();
	
	double minMeshSpacing=geometry->get_minMeshSpacing();

	double CFL_timeStep=CFL_number*minMeshSpacing/maxSoundSpeed;


	double maxDivergence=0.0;

	for (int i=1; i<=getSize().get_dir1(); ++i)
	for (int j=1; j<=getSize().get_dir2(); ++j)
	{
		double localDivergence=geometry->divergence(i,j, velocity); //units: 1/s
		if (fabs(localDivergence) > maxDivergence)
			maxDivergence=fabs(localDivergence);
	}

	double compression_timeStep=compression_number*(1.0/maxDivergence); //units: s
	

	return min_2arg(CFL_timeStep, compression_timeStep);
}

//***************************************************************************
//Calculate Data*************************************************************
double PerfectGas2D::calcPressure(int i, int j) const //units: Pa
{
	return molarDensity.get(i,j)*c_univGasConst*temperature.get(i,j);
		//units: Pa
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
