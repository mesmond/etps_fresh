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

class Operators2D : private StructuredGeometry2D
{
	private:
	double help;
	StructuredLocalField2D<double> area;
	StructuredLocalField2D<double> deltaFactor;
	StructuredLocalField2D<double> meshDelta;

	double localVolume;

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

	inline void print_geometry() const { this->print(); }
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
		//~ StructuredLocalField2D<double> cellFaceValues;
		//~ StructuredLocalField2D<double> localField;
	//~ 
		//~ Vector2D<double> gradient(0.0,0.0);
		//~ localField=scalarField.getLocalField(i,j);
		//~ deltaFactor=getLocalDeltaFactors(i,j);
//~ 
//~ 
		//~ cellFaceValues.N=localField.P*(1.0-deltaFactor.N)+deltaFactor.N*localField.N;
		//~ cellFaceValues.S=localField.P*(1.0-deltaFactor.S)+deltaFactor.S*localField.S;
		//~ cellFaceValues.E=localField.P*(1.0-deltaFactor.E)+deltaFactor.E*localField.E;
		//~ cellFaceValues.W=localField.P*(1.0-deltaFactor.W)+deltaFactor.W*localField.W;
//~ 
		//~ gradient.write_dir0((cellFaceValues.E-cellFaceValues.W)/deltaFactor.P_dir0);
		//~ gradient.write_dir1((cellFaceValues.N-cellFaceValues.S)/deltaFactor.P_dir1);
//~ 
		//~ return gradient;
	//~ }

	
};
