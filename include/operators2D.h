//***************************************************************************
//***************************************************************************
/*
 * File: operators2D.h
 * Description: This file contains class declarations for various
 * 	discretization operators in 2D.
 * 	The operators included are:
 * 		Gradient (cylindrical)
 * 		Divergence (cylindrical)
 * 		Gradient (cartesion)
 * 		Divergence (cartesion)
 *
 * 	The operators depend on the StructuredGeometry2D() class.
*/
//***************************************************************************
//***************************************************************************
#include <iostream>
#include <math.h>
#include <assert.h>

#include "dataTypes2D.h"
#include "mesmond-utils.h"

using namespace std;

class Operators2D : public StructuredGeometry2D
{
	private:
	double help;

	StructuredLocalField2D<double> deltaFactor;
	StructuredLocalField2D<double> meshDelta;

	double getVolume(int i, int j)
	{
		return 2.0*c_pi
			*this->getPoint(i,j).get_dir0()
			*this->getMeshDelta(i,j).get_dir0()
			*this->getMeshDelta(i,j).get_dir1();
	}

	StructuredLocalField2D<double> getCellAreas(int i, int j)
	{
		StructuredLocalField2D<double> area;

		area.N=2.0*c_pi
			*this->getPoint(i,j).get_dir0()
			*this->getMeshDelta(i,j).get_dir0();
		area.S=area.N;
		area.E=2.0*c_pi
			*(this->getPoint(i,j).get_dir0()
				+0.5*this->getMeshDelta(i,j).get_dir0())
			*this->getMeshDelta(i,j).get_dir1();
		area.W=2.0*c_pi
			*(this->getPoint(i,j).get_dir0()
				-0.5*this->getMeshDelta(i,j).get_dir0())
			*this->getMeshDelta(i,j).get_dir1();

		return area;
	}

	public:
	explicit Operators2D(
		Point2D<double> origin=Point2D<double>(0.0,0.0),
		Point2D<double> extent=Point2D<double>(1.0,1.0),
		Vector2D<int> size=Vector2D<int>(10,10),
		double refineMesh_dir0=1.0, double refineMesh_dir1=1.0)
		:	StructuredGeometry2D(origin, extent, size, refineMesh_dir0, refineMesh_dir1)
	{
		cout << "Constructor: Operators2D" << endl;
	}

	//~ inline void print_geometry() const { this->print(); }
	Vector2D<double> gradient(int i, int j, const SpacialArray2D<double>& scalarField)
	{
		this->checkIndices(i,j);

		StructuredLocalField2D<double> cellFaceValues;
		StructuredLocalField2D<double> localField;
	
		Vector2D<double> gradient(0.0,0.0);
		localField=scalarField.getLocalField(i,j);
		deltaFactor=getLocalDeltaFactors(i,j);


		cellFaceValues.N=localField.P*(1.0-deltaFactor.N)+deltaFactor.N*localField.N;
		cellFaceValues.S=localField.P*(1.0-deltaFactor.S)+deltaFactor.S*localField.S;
		cellFaceValues.E=localField.P*(1.0-deltaFactor.E)+deltaFactor.E*localField.E;
		cellFaceValues.W=localField.P*(1.0-deltaFactor.W)+deltaFactor.W*localField.W;

		gradient.write_dir0((cellFaceValues.E-cellFaceValues.W)/deltaFactor.P_dir0);
		gradient.write_dir1((cellFaceValues.N-cellFaceValues.S)/deltaFactor.P_dir1);

		return gradient;
	}

	//~ double divergence(int i, int j,
		//~ const SpacialArray2D<double>& scalar,
		//~ const SpacialArray2D<Vector2D<double> >& vector)
	//~ {
		//~ this->checkIndices(i,j);
//~ 
		//~ StructuredLocalField2D<double> scalarField;
		//~ StructuredLocalField2D<double> vectorField;
		//~ StructuredLocalField2D<double> localField;
		//~ StructuredLocalField2D<double> cellFaceValues;
	//~ 
		//~ scalarField=scalar.getLocalField(i,j);
		//~ vectorField=vector.getLocalField(i,j);
//~ 
		//~ localField.P_dir1=scalarField.P*vectorField.P.get_dir1();
		//~ localField.P_dir0=scalarField.P*vectorField.P.get_dir0();
		//~ 
		//~ localField.N=scalarField.N*vectorField.N.get_dir1();
		//~ localField.S=scalarField.S*vectorField.S.get_dir1();
		//~ localField.E=scalarField.E*vectorField.E.get_dir0();
		//~ localField.W=scalarField.W*vectorField.W.get_dir0();
		//~ 
		//~ deltaFactor=getLocalDeltaFactors(i,j);
//~ 
//~ 
		//~ cellFaceValues.N=localField.P_dir1*(1.0-deltaFactor.N)+deltaFactor.N*localField.N;
		//~ cellFaceValues.S=localField.P_dir1*(1.0-deltaFactor.S)+deltaFactor.S*localField.S;
		//~ cellFaceValues.E=localField.P_dir0*(1.0-deltaFactor.E)+deltaFactor.E*localField.E;
		//~ cellFaceValues.W=localField.P_dir0*(1.0-deltaFactor.W)+deltaFactor.W*localField.W;
//~ 
//~ 
//~ 
		//~ gradient.write_dir0((cellFaceValues.E-cellFaceValues.W)/deltaFactor.P_dir0);
		//~ gradient.write_dir1((cellFaceValues.N-cellFaceValues.S)/deltaFactor.P_dir1);
//~ 
		//~ return gradient;
	//~ }

	
};
