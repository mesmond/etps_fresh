//***************************************************************************
//***************************************************************************
/*
 * File: operators2D.cpp	
 * Description: This file contains the definitions for various
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
#include "operators2D.h"

Vector2D<double> Operators2D::gradient(int i, int j, const SpacialArray2D<double>& scalarField) const
{
	this->checkIndices(i,j);

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

	gradient.write_dir0((cellFaceValues.E-cellFaceValues.W)/deltaFactor.P_dir0);
	gradient.write_dir1((cellFaceValues.N-cellFaceValues.S)/deltaFactor.P_dir1);


	cout << "Gradient********************************" << endl;
	cout << "localField:" << endl;
	localField.print();
	
	cout << "cellFaces:" << endl;
	cellFaceValues.print();

	cout << "gradient=" << gradient << endl;
	

	return gradient;
}

double Operators2D::divergence(int i, int j,
	const SpacialArray2D<double>& scalar,
	const SpacialArray2D<Vector2D<double> >& vector) const
{
	this->checkIndices(i,j);

	StructuredLocalField2D<double> scalarField;
	StructuredLocalField2D<Vector2D<double> > vectorField;
	StructuredLocalField2D<double> localField;

	scalarField=scalar.getLocalField(i,j);
	vectorField=vector.getLocalField(i,j);

	localField.P_dir1=scalarField.P*vectorField.P.get_dir1();
	localField.P_dir0=scalarField.P*vectorField.P.get_dir0();
	
	localField.N=scalarField.N*vectorField.N.get_dir1();
	localField.S=scalarField.S*vectorField.S.get_dir1();
	localField.E=scalarField.E*vectorField.E.get_dir0();
	localField.W=scalarField.W*vectorField.W.get_dir0();



	cout << "Divergence (scalar*vector)****************" << endl;
	cout << "scalar:" << endl;
	scalarField.print();
	cout << "vector" << endl;
	vectorField.print();
	cout << "divergence=" << getDivergenceFromField(i,j,localField) << endl;




	
	return getDivergenceFromField(i,j,localField);
}

Vector2D<double> Operators2D::divergence(int i, int j,
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

	localField.P_dir0=dyadField.P.get_dir0();
	localField.P_dir1=dyadField.P.get_dir1();
	localField.N=dyadField.N.get_dir1();
	localField.S=dyadField.S.get_dir1();
	localField.E=dyadField.E.get_dir0();
	localField.W=dyadField.W.get_dir0();

	cout << "scalarField:" << endl;
	scalar_field.print();
	cout << "vec1_field:" << endl;
	vec1_field.print();
	cout << "vec2_field:" << endl;
	vec2_field.print();

	cout << "dyad Field:" << endl;
	dyadField.print();
	
	//~ localField.print();
	//~ cout << "localField.P=" << localField.P << endl; 

	cout << "divergence=" << getDivergenceFromField(i,j, localField) << endl;

	return getDivergenceFromField(i,j, localField);
}
