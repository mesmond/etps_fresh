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
#ifndef INCLUDE_OPERATORS2D_H_
#define INCLUDE_OPERATORS2D_H_

#include <iostream>
#include <math.h>
#include <assert.h>

#include "simulationConstants.h"
#include "dataTypes2D.h"

using namespace std;

class Operators2D : public StructuredGeometry2D
{
	private:
	double getVolume(int i, int j) const
	{
		return 2.0*c_pi
			*getPoint(i,j).get_dir1()
			*getMeshDelta(i,j).get_dir1()
			*getMeshDelta(i,j).get_dir2();
	}

	StructuredLocalField2D<double> getCellAreas(int i, int j) const
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

		return area;
	}

	template <typename T> T getDivergenceFromField(int i, int j,
		StructuredLocalField2D<T> field) const
	{
		StructuredLocalField2D<double> area;	
		StructuredLocalField2D<double> deltaFactor;
		StructuredLocalField2D<T> cellFlux;
		
		deltaFactor=getLocalDeltaFactors(i,j);
		area=getCellAreas(i,j);

		//Interpolate to cell faces.
		cellFlux.N=field.P_dir2*(1.0-deltaFactor.N)+deltaFactor.N*field.N;
		cellFlux.S=field.P_dir2*(1.0-deltaFactor.S)+deltaFactor.S*field.S;
		cellFlux.E=field.P_dir1*(1.0-deltaFactor.E)+deltaFactor.E*field.E;
		cellFlux.W=field.P_dir1*(1.0-deltaFactor.W)+deltaFactor.W*field.W;

		//Compute Divergence.
		T divergence=(1.0/getVolume(i,j))*(
			 cellNormal_N*cellFlux.N*area.N
			+cellNormal_S*cellFlux.S*area.S
			+cellNormal_E*cellFlux.E*area.E
			+cellNormal_W*cellFlux.W*area.W );

		return divergence;
	}

	public:
	explicit Operators2D(
		Point2D<double> origin=Point2D<double>(0.0,0.0),
		Point2D<double> extent=Point2D<double>(1.0,1.0),
		Vector2D<int> size=Vector2D<int>(10,10),
		double refineMesh_dir1=1.0, double refineMesh_dir2=1.0)
		:	StructuredGeometry2D(origin, extent, size, refineMesh_dir1, refineMesh_dir2)
	{
		cout << "Constructor: Operators2D" << endl;
	}

	Vector2D<double> gradient(const int i, const int j,
		const SpacialArray2D<double>& scalarField) const;
	double divergence(const int i, const int j,
		const SpacialArray2D<double>& scalar,
		const SpacialArray2D<Vector2D<double> >& vector) const;
	Vector2D<double> divergence(const int i, const int j,
		const SpacialArray2D<double>& scalar,
		const SpacialArray2D<Vector2D<double> >& vector1,
		const SpacialArray2D<Vector2D<double> >& vector2) const;
};


#endif // INCLUDE_OPERATORS2D_H_
