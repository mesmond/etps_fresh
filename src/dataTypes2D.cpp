//***************************************************************************
//***************************************************************************
/*
 * File: dataTypes2D.cpp
 * Description: This file contains class definitions for 2D array data types.
*/
//***************************************************************************
//***************************************************************************
#include "dataTypes2D.h"

//***************************************************************************
//***************************************************************************
//***************************************************************************
//Point2D********************************************************************

//***************************************************************************
//***************************************************************************
//***************************************************************************
//Vector2D*******************************************************************




//***************************************************************************
//***************************************************************************
//***************************************************************************
//Dyad2D*********************************************************************

//Define Operators***********************************************************
ostream& operator<<(ostream& os, const Dyad2D dyad)
{
	os << "[(" << dyad.value_00 << "," << dyad.value_01 << "),("
		<< dyad.value_10 << "," << dyad.value_11 << ")]";

	return os;
}

Dyad2D operator*(const Vector2D<double>& vec0, const Vector2D<double>& vec1)
{
	Dyad2D dyad;

	dyad.value_00=vec0.get_dir0()*vec1.get_dir0();
	dyad.value_01=vec0.get_dir0()*vec1.get_dir1();
	dyad.value_10=vec0.get_dir1()*vec1.get_dir0();
	dyad.value_11=vec0.get_dir1()*vec1.get_dir1();

	return dyad;
}



//***************************************************************************
//***************************************************************************
//***************************************************************************
//StructuredLocalField2D*****************************************************
template <typename type> void StructuredLocalField2D<type>::print() const
{
	cout << "StructuredLocalField2D:" << endl;
	cout << "\t P      =" << P << endl;
	cout << "\t P(dir0)=" << P_dir0 << endl;
	cout << "\t P(dir1)=" << P_dir1 << endl;
	cout << "\t N      =" << N << endl;
	cout << "\t S      =" << S << endl;
	cout << "\t E      =" << E << endl;
	cout << "\t W      =" << W << endl;
}

template class StructuredLocalField2D<Vector2D<double> >;
template class StructuredLocalField2D<double>;
template class StructuredLocalField2D<int>;

//***************************************************************************
//***************************************************************************
//***************************************************************************
//SpacialArray2D*************************************************************
template <typename type> SpacialArray2D<type>::SpacialArray2D( Vector2D<int> size )
{
	this->size=size;

	array = new type *[size.get_dir0()+2];
	for (int i=0; i<=size.get_dir0()+1; ++i)
	{
		array[i]= new type [size.get_dir1()+2];
	}

	for (int i=0; i<=size.get_dir0()+1; ++i) {
		for (int j=0; j<=size.get_dir1()+1; ++j)
		{
			array[i][j]=(type)0;
		}}
}

template <typename type> SpacialArray2D<type>::SpacialArray2D(const SpacialArray2D<type>& that)
{
	size=that.size;

	array = new type *[size.get_dir0()+2];
	for (int i=0; i<=size.get_dir0()+1; ++i)
	{
		array[i]= new type [size.get_dir1()+2];
	}

	for (int i=0; i<=size.get_dir0()+1; ++i) {
		for (int j=0; j<=size.get_dir1()+1; ++j)
		{
			array[i][j]=that.get(i,j);
		}}
}

template <typename type> SpacialArray2D<type>::~SpacialArray2D()
{
	for (int i=0; i<=size.get_dir0()+1; ++i)
	{
		delete[] array[i];
	}
	*array=0;
	delete array;
	array=0;
}


template <typename type> void SpacialArray2D<type>::print() const
{
	for (int j=size.get_dir1()+1; j>=0; --j)
	{
		for (int i=0; i<=size.get_dir0()+1; ++i)
		{
			cout << array[i][j] << " ";
		}
		cout << endl;
	}
}

template <typename type> void SpacialArray2D<type>::fill(const type& value)
{
	for (int i=0; i<=size.get_dir0()+1; ++i) {
	for (int j=0; j<=size.get_dir1()+1; ++j)
	{
		array[i][j]=value;
	}}
}


template <typename type> StructuredLocalField2D<type>
	SpacialArray2D<type>::getLocalField(int i, int j) const
{
	if (i < 1 || i > size.get_dir0() || j < 1 || j > size.get_dir1())
	{
		cout << "Array Index Out of range!" << endl;
		exit(1);

		return StructuredLocalField2D<type>();
	}
	else
	{
		StructuredLocalField2D<type> field;
		field.P=get(i,j);
		field.N=get(i,j+1);
		field.S=get(i,j-1);
		field.E=get(i+1,j);
		field.W=get(i-1,j);

		return field;
	}
}

template class SpacialArray2D<Vector2D<double> >;
template class SpacialArray2D<double>;
template class SpacialArray2D<int>;
//***************************************************************************
//***************************************************************************
//***************************************************************************
//StructuredGeometry2D*******************************************************
StructuredGeometry2D::StructuredGeometry2D(const StructuredGeometry2D& that)
{
	this->origin=that.origin;
	this->extent=that.extent;

	this->length=that.length;
	this->numZones=that.numZones;
	
	this->zoneDelta.dir0 = new double [numZones.get_dir0()+2];
	this->zoneDelta.dir1 = new double [numZones.get_dir1()+2];

	for (int i=0; i<=numZones.get_dir0()+1; ++i)
		this->zoneDelta.dir0[i]=that.zoneDelta.dir0[i];

	for (int j=0; j<=numZones.get_dir1()+1; ++j)
		this->zoneDelta.dir1[j]=that.zoneDelta.dir1[j];

	this->globalCoord.dir0 = new double [numZones.get_dir0()+2];
	this->globalCoord.dir1 = new double [numZones.get_dir1()+2];

	for (int i=0; i<=numZones.get_dir0()+1; ++i)
		this->globalCoord.dir0[i]=that.globalCoord.dir0[i];

	for (int j=0; j<=numZones.get_dir1()+1; ++j)
		this->globalCoord.dir1[j]=that.globalCoord.dir1[j];

	this->minMeshSpacing=that.minMeshSpacing;

}

StructuredGeometry2D::StructuredGeometry2D(
		Point2D<double> origin,
		Point2D<double> extent,
		Vector2D<int> size,
		double refineMesh_dir0,
		double refineMesh_dir1)
{
	this->origin=origin;
	this->extent=extent;

	length.write_dir0(extent.get_dir0()-origin.get_dir0());
	length.write_dir1(extent.get_dir1()-origin.get_dir1());

	assert(length.get_dir0() > 0.0);
	assert(length.get_dir1() > 0.0);

	numZones=size;

	zoneDelta.dir0 = new double [numZones.get_dir0()+2];// Add extra for ghost cells
	zoneDelta.dir1 = new double [numZones.get_dir1()+2];// Add extra for ghost cells

	globalCoord.dir0 = new double [numZones.get_dir0()+2];// Add extra for ghost cells
	globalCoord.dir1 = new double [numZones.get_dir1()+2];// Add extra for ghost cells

	//Mesh Spacing***********************************************************
	//Given a mesh refinement constant, each cell in the specified direction
	//	will be a factor of the previous cell in that direction.
	//	To establish the mesh spacings, we first solve for the origin mesh
	//	size.
	//Get origin mesh spacing.
	double sum_dir0=0.0;
	for (int i=0; i<numZones.get_dir0(); ++i)
		sum_dir0+=pow(refineMesh_dir0, (double)i);

	//Set mesh spacing values
	//The ghost cells dimensions are set to the
	//	boundary cell dimenstions.  These must be shared with
	//	neighboring geometries when applicable.
	double initialSpacing_dir0=length.get_dir0()/sum_dir0;
	zoneDelta.dir0[0]=initialSpacing_dir0;
	zoneDelta.dir0[1]=zoneDelta.dir0[0];
	for (int i=2; i<=numZones.get_dir0(); ++i)
		zoneDelta.dir0[i]=pow(refineMesh_dir0,i-1)*zoneDelta.dir0[1];
	zoneDelta.dir0[numZones.get_dir0()+1]=zoneDelta.dir0[numZones.get_dir0()];

	//Get origin mesh spacing.
	double sum_dir1=0.0;
	for (int i=0; i<numZones.get_dir1(); ++i)
		sum_dir1+=pow(refineMesh_dir1, (double)i);

	//Set mesh spacing values
	//The ghost cells dimensions are set to the
	//	boundary cell dimenstions.  These must be shared with
	//	neighboring geometries when applicable.
	double initialSpacing_dir1=length.get_dir1()/sum_dir1;
	zoneDelta.dir1[0]=initialSpacing_dir1;
	zoneDelta.dir1[1]=zoneDelta.dir1[0];
	for (int j=2; j<=numZones.get_dir1(); ++j)
		zoneDelta.dir1[j]=pow(refineMesh_dir1,j-1)*zoneDelta.dir1[1];
	zoneDelta.dir1[numZones.get_dir1()+1]=zoneDelta.dir1[numZones.get_dir1()];

	//Global Cell Coordinates************************************************
	double location=origin.get_dir0();
	globalCoord.dir0[0]=location-0.5*zoneDelta.dir0[0];
	for (int i=1; i<=numZones.get_dir0()+1; ++i)
	{
		location+=zoneDelta.dir0[i];
		globalCoord.dir0[i]=location-0.5*zoneDelta.dir0[i];
	}

	location=origin.get_dir1();
	globalCoord.dir1[0]=location-0.5*zoneDelta.dir1[0];
	for (int j=1; j<=numZones.get_dir1()+1; ++j)
	{
		location+=zoneDelta.dir1[j];
		globalCoord.dir1[j]=location-0.5*zoneDelta.dir1[j];
	}

	//Define minMeshSpacing**********************************************
	double minValue=length.get_dir0();
	for (int i=1; i<=numZones.get_dir0(); ++i)
		if (zoneDelta.dir0[i] < minValue) minValue=zoneDelta.dir0[i];

	for (int j=1; j<=numZones.get_dir1(); ++j)
		if (zoneDelta.dir1[j] < minValue) minValue=zoneDelta.dir1[j];

	minMeshSpacing=minValue;
}

StructuredGeometry2D::~StructuredGeometry2D()
{
	delete[] zoneDelta.dir0;
	delete[] zoneDelta.dir1;

	zoneDelta.dir0=NULL;
	zoneDelta.dir1=NULL;

	delete[] globalCoord.dir0;
	delete[] globalCoord.dir1;

	globalCoord.dir0=NULL;
	globalCoord.dir1=NULL;
}

void StructuredGeometry2D::print() const
{
	int maxi=numZones.get_dir0()+1;
	int maxj=numZones.get_dir1()+1;
	
	cout << "StructuredGeometry2D: \n"
		<< "\t Location: " << origin << "->" << extent << endl
		<< "\t Size=" << numZones << endl;
	cout << "\t Cell Center Locations:" << endl;

	for (int j=maxj; j>=0; --j)
	{
		cout << "\t" << globalCoord.dir1[j] << endl;
	}

	cout << "\t\t";
	for (int i=0; i<=maxi; ++i)
	{
		cout << globalCoord.dir0[i] << " ";
	}
	cout << endl;

}

void StructuredGeometry2D::link_north(const StructuredGeometry2D& toLink)
{
	assert(getCount_dir0() == toLink.getCount_dir0());
	assert( isEqual(get_extent().get_dir1(), toLink.get_origin().get_dir1()));
	setBdy_zoneDelta_north(toLink.get_zoneDelta_south());
}

void StructuredGeometry2D::link_south(const StructuredGeometry2D& toLink)
{
	assert(getCount_dir0() == toLink.getCount_dir0());
	assert( isEqual(get_origin().get_dir1(), toLink.get_extent().get_dir1()));
	setBdy_zoneDelta_south(toLink.get_zoneDelta_north());
}

void StructuredGeometry2D::link_east(const StructuredGeometry2D& toLink)
{
	assert(getCount_dir1() == toLink.getCount_dir1());
	assert( isEqual(get_extent().get_dir0(), toLink.get_origin().get_dir0()));
	setBdy_zoneDelta_east(toLink.get_zoneDelta_west());
}

void StructuredGeometry2D::link_west(const StructuredGeometry2D& toLink)
{
	assert(getCount_dir1() == toLink.getCount_dir1());
	assert( isEqual(get_origin().get_dir0(), toLink.get_extent().get_dir0()));
	setBdy_zoneDelta_west(toLink.get_zoneDelta_east());
}



//***************************************************************************
//***************************************************************************
//***************************************************************************
//spacialArray2D*************************************************************
//~ spacialArray2D::spacialArray2D(const spacialArray2D& that)
//~ {
	//~ this->numZones.dir0=that.numZones.dir0;
	//~ this->numZones.dir1=that.numZones.dir1;
	//~ 
	//~ array= new double *[numZones.dir0+2];
	//~ for (int i=0; i<=numZones.dir0+1; ++i)
	//~ {
		//~ array[i]= new double [numZones.dir1+2];
	//~ }
//~ 
	//~ for (int i=0; i<=numZones.dir0+1; ++i) {
	//~ for (int j=0; j<=numZones.dir1+1; ++j)
	//~ {
		//~ array[i][j]= that.array[i][j];
	//~ }}
//~ }
//~ 
//~ spacialArray2D::spacialArray2D(int size_dir0, int size_dir1, double valueInit)
//~ {
	//~ numZones.dir0=size_dir0;
	//~ numZones.dir1=size_dir1;
//~ 
	//~ array= new double *[numZones.dir0+2];
	//~ for (int i=0; i<=numZones.dir0+1; ++i)
	//~ {
		//~ array[i]= new double [numZones.dir1+2];
	//~ }
//~ 
	//~ for (int i=0; i<=numZones.dir0+1; ++i) {
	//~ for (int j=0; j<=numZones.dir1+1; ++j)
	//~ {
		//~ array[i][j] = valueInit;
	//~ }}
//~ 
	//~ //cout << "Done" << endl;
//~ }
//~ 
//~ spacialArray2D& spacialArray2D::operator=(const spacialArray2D& that)
//~ {
	//~ this->numZones.dir0=that.numZones.dir0;
	//~ this->numZones.dir1=that.numZones.dir1;
	//~ 
	//~ array= new double *[numZones.dir0+2];
	//~ for (int i=0; i<=numZones.dir0+1; ++i)
	//~ {
		//~ array[i]= new double [numZones.dir1+2];
	//~ }
//~ 
	//~ for (int i=0; i<=numZones.dir0+1; ++i) {
	//~ for (int j=0; j<=numZones.dir1+1; ++j)
	//~ {
		//~ array[i][j]= that.array[i][j];
	//~ }}
//~ 
	//~ return *this;
//~ }
//~ 
//~ spacialArray2D::~spacialArray2D()
//~ {
	//~ for (int i=0; i<=numZones.dir0+1; ++i)
	//~ {
		//~ delete [] array[i];
//~ 
		//~ *array=NULL;
	//~ }
	//~ delete [] array;
//~ 
	//~ array=NULL;
//~ }
//~ 
//~ double spacialArray2D::get(int i, int j) const
//~ {
	//~ return array[i][j];
//~ }
//~ void spacialArray2D::write(int i, int j, double value)
//~ {
	//~ array[i][j]=value;
//~ }
//~ int spacialArray2D::getSize_dir0() const
//~ {
	//~ return numZones.dir0;
//~ }
//~ int spacialArray2D::getSize_dir1() const	
//~ {
	//~ return numZones.dir1;
//~ }

//~ cellDirections spacialArray2D::getLocalField(int i, int j) const
//~ {
	//~ cellDirections field;
	//~ field.P=array[i][j];
	//~ field.N=array[i][j+1];
	//~ field.S=array[i][j-1];
	//~ field.E=array[i+1][j];
	//~ field.W=array[i-1][j];
//~ 
	//~ return field;
//~ }

//***************************************************************************
//***************************************************************************
//***************************************************************************
//vectorArray2D**************************************************************
//~ double vectorArray2D::getMagnitude(int i, int j)
//~ {
	//~ return sqrt(
		//~ pow(array_dir0.get(i,j),2.0)
		//~ +pow(array_dir1.get(i,j),2.0));
//~ }
//~ 
//~ void vectorArray2D::writeComponent_dir0(int i, int j, double value)
//~ {
	//~ array_dir0.write(i,j,value);
//~ }
//~ void vectorArray2D::writeComponent_dir1(int i, int j, double value)
//~ {
	//~ array_dir1.write(i,j,value);
//~ }
//~ 
//~ double vectorArray2D::component_dir0(int i, int j)
//~ {
	//~ return array_dir0.get(i,j);
//~ }
//~ double vectorArray2D::component_dir1(int i, int j)
//~ {
	//~ return array_dir1.get(i,j);
//~ }
//~ 
//~ int vectorArray2D::getSize_dir0()
//~ {
	//~ return array_dir1.getSize_dir0();
//~ }
//~ int vectorArray2D::getSize_dir1()
//~ {
	//~ return array_dir1.getSize_dir1();
//~ }

/*
 * Function: vectorArray2D::product()
 * Purpose: Compute and return a local field determined from the
 * 	product of a scalar and a vector.
 *
 * Description: The product of a scalar and a vector yield a vector
 * 	flux.  A local vector flux field is computed here and returned.
 */
//~ cellDirections vectorArray2D::product(
	//~ cellDirections scalarField,
	//~ int i, int j) const
//~ {
	//~ cellDirections fieldFlux;
	//~ fieldFlux.P_dir0=scalarField.P*array_dir0.get(i,j);
	//~ fieldFlux.P_dir1=scalarField.P*array_dir1.get(i,j);
//~ 
	//~ fieldFlux.N=scalarField.N*array_dir1.get(i,j+1);
	//~ fieldFlux.S=scalarField.S*array_dir1.get(i,j-1);
	//~ fieldFlux.E=scalarField.E*array_dir0.get(i+1,j);
	//~ fieldFlux.W=scalarField.W*array_dir0.get(i-1,j);
//~ 
	//~ return fieldFlux;
//~ }

