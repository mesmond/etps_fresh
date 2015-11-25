//***************************************************************************
//***************************************************************************
/*
 * File: structuredGeometry2D.cpp
 * Description: This file contains definitions for the
 * 	StructuredGeometry2D class.  It contains not only
 * 	mesh information, setters and getters, but also differential
 * 	operators such as
 * 		gradient (scalar)
 * 		divergence (scalar, vector)
 * 		divergence (scalar, vector, vector)
*/
//***************************************************************************
//***************************************************************************
#include "structuredGeometry2D.h"

//***************************************************************************
//***************************************************************************
//***************************************************************************
//StructuredGeometry2D*******************************************************
StructuredGeometry2D::StructuredGeometry2D(
		Point2D<double> origin,
		Point2D<double> extent,
		Vector2D<int> size,
		Vector2D<double> refineMesh)
{
	this->origin=origin;
	this->extent=extent;

	length.write_dir1(extent.get_dir1()-origin.get_dir1());
	length.write_dir2(extent.get_dir2()-origin.get_dir2());

	this->refineMesh=refineMesh;

	assert(length.get_dir1() > 0.0);
	assert(length.get_dir2() > 0.0);

	numZones=size;

	zoneDelta.dir1 = new double [numZones.get_dir1()+2];// Add extra for ghost cells
	zoneDelta.dir2 = new double [numZones.get_dir2()+2];// Add extra for ghost cells

	globalCoord.dir1 = new double [numZones.get_dir1()+2];// Add extra for ghost cells
	globalCoord.dir2 = new double [numZones.get_dir2()+2];// Add extra for ghost cells

	//Mesh Spacing***********************************************************
	//Given a mesh refinement constant, each cell in the specified direction
	//	will be a factor of the previous cell in that direction.
	//	To establish the mesh spacings, we first solve for the origin mesh
	//	size.
	//Get origin mesh spacing.
	double sum_dir1=0.0;
	for (int i=0; i<numZones.get_dir1(); ++i)
		sum_dir1+=pow(refineMesh.get_dir1(), (double)i);

	//Set mesh spacing values
	//The ghost cells dimensions are set to the
	//	boundary cell dimenstions.  These must be shared with
	//	neighboring geometries when applicable.
	double initialSpacing_dir1=length.get_dir1()/sum_dir1;
	zoneDelta.dir1[0]=initialSpacing_dir1;
	zoneDelta.dir1[1]=zoneDelta.dir1[0];
	for (int i=2; i<=numZones.get_dir1(); ++i)
		zoneDelta.dir1[i]=pow(refineMesh.get_dir1(),i-1)*zoneDelta.dir1[1];
	zoneDelta.dir1[numZones.get_dir1()+1]=zoneDelta.dir1[numZones.get_dir1()];

	//Get origin mesh spacing.
	double sum_dir2=0.0;
	for (int i=0; i<numZones.get_dir2(); ++i)
		sum_dir2+=pow(refineMesh.get_dir2(), (double)i);

	//Set mesh spacing values
	//The ghost cells dimensions are set to the
	//	boundary cell dimenstions.  These must be shared with
	//	neighboring geometries when applicable.
	double initialSpacing_dir2=length.get_dir2()/sum_dir2;
	zoneDelta.dir2[0]=initialSpacing_dir2;
	zoneDelta.dir2[1]=zoneDelta.dir2[0];
	for (int j=2; j<=numZones.get_dir2(); ++j)
		zoneDelta.dir2[j]=pow(refineMesh.get_dir2(),j-1)*zoneDelta.dir2[1];
	zoneDelta.dir2[numZones.get_dir2()+1]=zoneDelta.dir2[numZones.get_dir2()];

	//Global Cell Coordinates************************************************
	double location=origin.get_dir1();
	globalCoord.dir1[0]=location-0.5*zoneDelta.dir1[0];
	for (int i=1; i<=numZones.get_dir1()+1; ++i)
	{
		location+=zoneDelta.dir1[i];
		globalCoord.dir1[i]=location-0.5*zoneDelta.dir1[i];
	}

	location=origin.get_dir2();
	globalCoord.dir2[0]=location-0.5*zoneDelta.dir2[0];
	for (int j=1; j<=numZones.get_dir2()+1; ++j)
	{
		location+=zoneDelta.dir2[j];
		globalCoord.dir2[j]=location-0.5*zoneDelta.dir2[j];
	}

	//Define minMeshSpacing**********************************************
	double minValue=length.get_dir1();
	for (int i=1; i<=numZones.get_dir1(); ++i)
		if (zoneDelta.dir1[i] < minValue) minValue=zoneDelta.dir1[i];

	for (int j=1; j<=numZones.get_dir2(); ++j)
		if (zoneDelta.dir2[j] < minValue) minValue=zoneDelta.dir2[j];

	minMeshSpacing=minValue;
}

StructuredGeometry2D::StructuredGeometry2D(const StructuredGeometry2D& that)
{
	this->origin=that.origin;
	this->extent=that.extent;

	this->length=that.length;
	this->numZones=that.numZones;
	
	this->zoneDelta.dir1 = new double [numZones.get_dir1()+2];
	this->zoneDelta.dir2 = new double [numZones.get_dir2()+2];

	for (int i=0; i<=numZones.get_dir1()+1; ++i)
		this->zoneDelta.dir1[i]=that.zoneDelta.dir1[i];

	for (int j=0; j<=numZones.get_dir2()+1; ++j)
		this->zoneDelta.dir2[j]=that.zoneDelta.dir2[j];

	this->globalCoord.dir1 = new double [numZones.get_dir1()+2];
	this->globalCoord.dir2 = new double [numZones.get_dir2()+2];

	for (int i=0; i<=numZones.get_dir1()+1; ++i)
		this->globalCoord.dir1[i]=that.globalCoord.dir1[i];

	for (int j=0; j<=numZones.get_dir2()+1; ++j)
		this->globalCoord.dir2[j]=that.globalCoord.dir2[j];

	this->minMeshSpacing=that.minMeshSpacing;

}

StructuredGeometry2D::~StructuredGeometry2D()
{
	delete[] zoneDelta.dir1;
	delete[] zoneDelta.dir2;

	zoneDelta.dir1=nullptr;
	zoneDelta.dir2=nullptr;

	delete[] globalCoord.dir1;
	delete[] globalCoord.dir2;

	globalCoord.dir1=nullptr;
	globalCoord.dir2=nullptr;
}

//Output*********************************************************************
void StructuredGeometry2D::print() const
{
	int maxi=numZones.get_dir1()+1;
	int maxj=numZones.get_dir2()+1;
	
	cout << "StructuredGeometry2D: \n"
		<< "\t Location: " << origin << "->" << extent << endl
		<< "\t Size=" << numZones << endl;
	cout << "\t Cell Center Locations:" << endl;

	for (int j=maxj; j>=0; --j)
	{
		cout << "\t" << globalCoord.dir2[j] << endl;
	}

	cout << "\t\t";
	for (int i=0; i<=maxi; ++i)
	{
		cout << globalCoord.dir1[i] << " ";
	}
	cout << endl;

}


void StructuredGeometry2D::vtkOutput(const char* fileName) const
{
	cout << "Preparing Output to VTK" << endl;

	//Store the output precision so that it can be returned to its
	//	usual value after this function.
	streamsize ss = cout.precision();


	//stringstream outputIndex;
	//outputIndex << setw(5) << setfill('0') << outputCount;
	//string numberString = outputIndex.str();

	//~ string outFileName=inputFilePrefix+"_block#"+to_string(blockID)+"@-"
		//~ +numberString+".vtk";

	ofstream output;
	output.open(fileName, ios::out);


	if ( output.is_open() )
	{
		int maxi=numZones.get_dir1()+1;
		int maxj=numZones.get_dir2()+1;
		output << "# vtk DataFile Version 2.0" << endl;
		output 	<< "Blank"
				<< scientific
				<< setprecision(10)
				<< endl;
		output << "ASCII" << endl;
		output	<< "DATASET STRUCTURED_GRID" << endl;

		output << "DIMENSIONS " << maxi << " " 
				<< maxj << " 2" << endl;

		output << "POINTS " <<
			maxi*maxj*2
			<< " float\n";

		for (double depth=0.0; depth<0.00101; depth=depth+0.001)
		{
			for (int j=1; j<=maxj; ++j) {
			for (int i=1; i<=maxi; ++i)
			{
				output	<< globalCoord.dir1[i]-0.5*zoneDelta.dir1[i]
						<< " "
						<< globalCoord.dir2[j]-0.5*zoneDelta.dir2[j]
						<< " "
						<< depth << endl;
			}}
		}

		output.close();

	}

	//Restore original cout precision.
	cout.precision (ss);
	cout.unsetf (ios::floatfield); // unsets fixed
}

void StructuredGeometry2D::vtkOutput(const SpacialArray2D<double>& scalar,
	const char* fileName) const
{
	vtkOutput(fileName);

	//Store the output precision so that it can be returned to its
	//	usual value after this function.
	streamsize ss = cout.precision();

	assert(scalar.getCount_dir1() == numZones.get_dir1());
	assert(scalar.getCount_dir2() == numZones.get_dir2());

	ofstream output;
	output.open(fileName, ios::app);

	if ( output.is_open() )
	{
		output << "\nCELL_DATA "
			<< numZones.get_dir2()*numZones.get_dir1() << endl;
		output << "SCALARS " << "scalar" << " float" << endl;
		output << "LOOKUP_TABLE default" << endl;
		for (int j=1; j<=numZones.get_dir2(); ++j){
		for (int i=1; i<=numZones.get_dir1(); ++i)
		{
			output << scalar.get(i,j) << endl;
		}}
		output << endl;
		
		output.close();
	}

	//Restore original cout precision.
	cout.precision (ss);
	cout.unsetf (ios::floatfield); // unsets fixed
}

void StructuredGeometry2D::vtkOutput(
	const SpacialArray2D<double>& scalar,
	const SpacialArray2D<Vector2D<double> >& vector,
	const char* fileName) const
{
	vtkOutput(scalar, fileName);

	//Store the output precision so that it can be returned to its
	//	usual value after this function.
	streamsize ss = cout.precision();

	assert(vector.getCount_dir1() == numZones.get_dir1());
	assert(vector.getCount_dir2() == numZones.get_dir2());

	ofstream output;
	output.open(fileName, ios::app);

	if ( output.is_open() )
	{
		output << "VECTORS " << "vector" << " float" << endl;
		for (int j=1; j<=numZones.get_dir2(); ++j){
		for (int i=1; i<=numZones.get_dir1(); ++i)
		{
			output << vector.get(i,j).get_dir1()
				<< " " << vector.get(i,j).get_dir2()
				<< " 0.0"
				<< endl;
		}}
		output << endl;
		
		output.close();
	}

	//Restore original cout precision.
	cout.precision (ss);
	cout.unsetf (ios::floatfield); // unsets fixed
}

void StructuredGeometry2D::vtkOutput(
	const SpacialArray2D<double>& scalar,
	const SpacialArray2D<Vector2D<double> >& vector,
	const SpacialArray2D<double>& divergence,
	const char* fileName) const
{
	vtkOutput(scalar, vector, fileName);

	//Store the output precision so that it can be returned to its
	//	usual value after this function.
	streamsize ss = cout.precision();

	assert(divergence.getCount_dir1() == numZones.get_dir1());
	assert(divergence.getCount_dir2() == numZones.get_dir2());

	ofstream output;
	output.open(fileName, ios::app);

	if ( output.is_open() )
	{
		output << "\nSCALARS " << "divergence" << " float" << endl;
		output << "LOOKUP_TABLE default" << endl;
		for (int j=1; j<=numZones.get_dir2(); ++j){
		for (int i=1; i<=numZones.get_dir1(); ++i)
		{
			output << divergence.get(i,j) << endl;
		}}
		output << endl;
		
		output.close();
	}

	//Restore original cout precision.
	cout.precision (ss);
	cout.unsetf (ios::floatfield); // unsets fixed
}

//Diferential Operators******************************************************
Vector2D<double> StructuredGeometry2D::gradient(int i, int j,
	const SpacialArray2D<double>& scalarField) const
{
	checkIndices(i,j);

	StructuredLocalField2D<double> cellFaceValues;
	StructuredLocalField2D<double> localField;
	StructuredLocalField2D<double> deltaFactor;

	Vector2D<double> gradient(0.0,0.0);
	localField=scalarField.getLocalField(i,j);
	deltaFactor=getLocalDeltaFactors(i,j);


	cellFaceValues.N=localField.P*(1.0-deltaFactor.N)+deltaFactor.N*localField.N;
	cellFaceValues.S=localField.P*(1.0-deltaFactor.S)+deltaFactor.S*localField.S;
	cellFaceValues.E=localField.P*(1.0-deltaFactor.E)+deltaFactor.E*localField.E;
	cellFaceValues.W=localField.P*(1.0-deltaFactor.W)+deltaFactor.W*localField.W;

	gradient.write_dir1((cellFaceValues.E-cellFaceValues.W)/deltaFactor.P_dir1);
	gradient.write_dir2((cellFaceValues.N-cellFaceValues.S)/deltaFactor.P_dir2);

	return gradient;
}


double StructuredGeometry2D::divergence(int i, int j,
	const SpacialArray2D<double>& scalar,
	const SpacialArray2D<Vector2D<double> >& vector) const
{
	this->checkIndices(i,j);

	StructuredLocalField2D<double> scalarField;
	StructuredLocalField2D<Vector2D<double> > vectorField;
	StructuredLocalField2D<double> localField;

	scalarField=scalar.getLocalField(i,j);
	vectorField=vector.getLocalField(i,j);

	localField.P_dir2=scalarField.P*vectorField.P.get_dir2();
	localField.P_dir1=scalarField.P*vectorField.P.get_dir1();
	
	localField.N=scalarField.N*vectorField.N.get_dir2();
	localField.S=scalarField.S*vectorField.S.get_dir2();
	localField.E=scalarField.E*vectorField.E.get_dir1();
	localField.W=scalarField.W*vectorField.W.get_dir1();
	
	return getDivergenceFromField(i,j,localField);
}


Vector2D<double> StructuredGeometry2D::divergence(int i, int j,
	const SpacialArray2D<double>& scalar,
	const SpacialArray2D<Vector2D<double> >& vector1,
	const SpacialArray2D<Vector2D<double> >& vector2) const
{
	this->checkIndices(i,j);
	Vector2D<double> result(0.0,0.0);

	StructuredLocalField2D<double> scalar_field;
	StructuredLocalField2D<Vector2D<double> > vec1_field;
	StructuredLocalField2D<Vector2D<double> > vec2_field;
	StructuredLocalField2D<Dyad2D> dyadField;
	StructuredLocalField2D<Vector2D<double> > localField;

	scalar_field=scalar.getLocalField(i,j);
	vec1_field=vector1.getLocalField(i,j);
	vec2_field=vector2.getLocalField(i,j);

	dyadField.P=scalar_field.P*(vec1_field.P*vec2_field.P);
	dyadField.N=scalar_field.N*(vec1_field.N*vec2_field.N);
	dyadField.S=scalar_field.S*(vec1_field.S*vec2_field.S);
	dyadField.E=scalar_field.E*(vec1_field.E*vec2_field.E);
	dyadField.W=scalar_field.W*(vec1_field.W*vec2_field.W);

	localField.P_dir1=dyadField.P.get_dir1();
	localField.P_dir2=dyadField.P.get_dir2();
	localField.N=dyadField.N.get_dir2();
	localField.S=dyadField.S.get_dir2();
	localField.E=dyadField.E.get_dir1();
	localField.W=dyadField.W.get_dir1();

	return getDivergenceFromField(i,j, localField);
}








//Mutation Handling**********************************************************
void StructuredGeometry2D::link_north(const StructuredGeometry2D& toLink)
{
	assert(getCount_dir1() == toLink.getCount_dir1());
	assert( isEqual(get_extent().get_dir2(), toLink.get_origin().get_dir2()));
	setBdy_zoneDelta_north(toLink.get_zoneDelta_south());
}

void StructuredGeometry2D::link_south(const StructuredGeometry2D& toLink)
{
	assert(getCount_dir1() == toLink.getCount_dir1());
	assert( isEqual(get_origin().get_dir2(), toLink.get_extent().get_dir2()));
	setBdy_zoneDelta_south(toLink.get_zoneDelta_north());
}

void StructuredGeometry2D::link_east(const StructuredGeometry2D& toLink)
{
	assert(getCount_dir2() == toLink.getCount_dir2());
	assert( isEqual(get_extent().get_dir1(), toLink.get_origin().get_dir1()));
	setBdy_zoneDelta_east(toLink.get_zoneDelta_west());
}

void StructuredGeometry2D::link_west(const StructuredGeometry2D& toLink)
{
	assert(getCount_dir2() == toLink.getCount_dir2());
	assert( isEqual(get_origin().get_dir1(), toLink.get_extent().get_dir1()));
	setBdy_zoneDelta_west(toLink.get_zoneDelta_east());
}

//Get Information************************************************************
StructuredLocalField2D<double>
	StructuredGeometry2D::getLocalDeltaFactors(int i, int j) const
{
	StructuredLocalField2D<double> result;

	result.N=zoneDelta.dir2[j]/(zoneDelta.dir2[j]+zoneDelta.dir2[j+1]);
	result.S=zoneDelta.dir2[j]/(zoneDelta.dir2[j]+zoneDelta.dir2[j-1]);
	result.E=zoneDelta.dir1[i]/(zoneDelta.dir1[i]+zoneDelta.dir1[i+1]);
	result.W=zoneDelta.dir1[i]/(zoneDelta.dir1[i]+zoneDelta.dir1[i-1]);

	result.P_dir1=zoneDelta.dir1[i];
	result.P_dir2=zoneDelta.dir2[j];

	return result;
}

Point2D<double> StructuredGeometry2D::getPoint(int i, int j) const
{
	return Point2D<double>(globalCoord.dir1[i], globalCoord.dir2[j]);
}

Vector2D<double> StructuredGeometry2D::getMeshDelta(int i, int j) const
{
	return Vector2D<double>(zoneDelta.dir1[i], zoneDelta.dir2[j]);
}

//Coordinate System Specific Functions: Default to rectangular.
double StructuredGeometry2D::getVolume(int i, int j) const
{
	double unitWidth=1.0;
	return unitWidth
		*getMeshDelta(i,j).get_dir1()*getMeshDelta(i,j).get_dir2();
}

StructuredLocalField2D<double> StructuredGeometry2D::getCellAreas(int i, int j) const
{
	double unitWidth=1.0;
	StructuredLocalField2D<double> area;
	area.N=unitWidth*getMeshDelta(i,j).get_dir1();
	area.S=area.N;
	area.E=unitWidth*getMeshDelta(i,j).get_dir2();
	area.W=area.E;

	area.P=0.0;
	area.P_dir1=0.0;
	area.P_dir2=0.0;

	return area;
}




//Set Information********************************************************
void StructuredGeometry2D::setBdy_zoneDelta_north(double value)
{
	int index=numZones.get_dir2()+1;
	zoneDelta.dir2[index] = value;
	globalCoord.dir2[index]=extent.get_dir2()+0.5*zoneDelta.dir2[index];
}

void StructuredGeometry2D::setBdy_zoneDelta_south(double value)
{
	int index=0;
	zoneDelta.dir2[index] = value;
	globalCoord.dir2[index]=origin.get_dir2()-0.5*zoneDelta.dir2[index];
}

void StructuredGeometry2D::setBdy_zoneDelta_east(double value)
{
	int index=numZones.get_dir1()+1;
	zoneDelta.dir1[index] = value;
	globalCoord.dir1[index]=extent.get_dir1()+0.5*zoneDelta.dir1[index];
}

void StructuredGeometry2D::setBdy_zoneDelta_west(double value)
{
	int index=0;
	zoneDelta.dir1[index] = value;
	globalCoord.dir1[index]=origin.get_dir1()-0.5*zoneDelta.dir1[index];
}

//Check Information**********************************************************
void StructuredGeometry2D::checkIndices(int i, int j) const
{
	assert( i > 0 );
	assert( j > 0 );
	
	assert( i <= numZones.get_dir1() );
	assert( j <= numZones.get_dir2() );
}



//****************************************************************Cylindrical
//****************************************************************Cylindrical
//****************************************************************Cylindrical
//StructuredCylGeometry2D*****************************************Cylindrical
double StructuredCylGeometry2D::getVolume(int i, int j) const
{
	return 2.0*c_pi
		*getPoint(i,j).get_dir1()
		*getMeshDelta(i,j).get_dir1()
		*getMeshDelta(i,j).get_dir2();
}

StructuredLocalField2D<double> StructuredCylGeometry2D::getCellAreas(int i, int j) const
{
	StructuredLocalField2D<double> area;

	area.N=2.0*c_pi
		*getPoint(i,j).get_dir1()
		*getMeshDelta(i,j).get_dir1();
	area.S=area.N;
	area.E=2.0*c_pi
		*(getPoint(i,j).get_dir1()
			+0.5*getMeshDelta(i,j).get_dir1())
		*getMeshDelta(i,j).get_dir2();
	area.W=2.0*c_pi
		*(getPoint(i,j).get_dir1()
			-0.5*getMeshDelta(i,j).get_dir1())
		*getMeshDelta(i,j).get_dir2();

	area.P=0.0;
	area.P_dir1=0.0;
	area.P_dir2=0.0;

	return area;
}





