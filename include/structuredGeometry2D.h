//***************************************************************************
//***************************************************************************
/*
 * File: structuredGeometry2D.h
 * Description: This file contains class declarations for the
 * 	StructuredGeometry2D class.  It contains not only
 * 	mesh information, setters and getters, but also differential
 * 	operators such as
 * 		gradient (scalar)
 * 		divergence (scalar, vector)
 * 		divergence (scalar, vector, vector)
*/
//***************************************************************************
//***************************************************************************
#ifndef INCLUDE_STRUCTUREDGEOMETRY2D_H_
#define INCLUDE_STRUCTUREDGEOMETRY2D_H_

#include <iostream>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <iomanip>

#include "simulationConstants.h"
#include "mesmond-utils.h"
#include "dataTypes2D.h"

using namespace std;

//***************************************************************************
//***************************************************************************
//***************************************************************************
//StructuredGeometry2D*******************************************************
/*
 * class: StructuredGeometry2D
 * Purpose: Hold and process geometric parameters for a structured geometry.
 *
 * Description: StructuredGeometry2D holds information on the spacial placement
 * 	of a particular geometry. It also holds information
 * 	regarding the mesh spacing and mesh density.
 * 	Global Coordinates of structured cells
 * 	are specified using a cell-centered scheme.
 * 	This class is designed for an orthogonal grid.
 */
class StructuredGeometry2D
{
	private:
	Point2D<double> origin;
	Point2D<double> extent;
	Vector2D<double> length;

	Vector2D<double> refineMesh;

	Vector2D<int> numZones;

	struct MeshSpacing
	{
		double *dir1;
		double *dir2;
	} zoneDelta, globalCoord;

	double minMeshSpacing;

	public:
	explicit StructuredGeometry2D(
		Point2D<double> origin=Point2D<double>(0.0,0.0),
		Point2D<double> extent=Point2D<double>(1.0,1.0),
		Vector2D<int> size=Vector2D<int>(10,10),
		Vector2D<double> refineMesh=Vector2D<double>(1.0,1.0)); // Con
	StructuredGeometry2D(const StructuredGeometry2D& that); // Copy Con
	StructuredGeometry2D(StructuredGeometry2D&& other) = delete; //Move Con
	StructuredGeometry2D& operator=(StructuredGeometry2D&& other) = delete; //Move Assign
	StructuredGeometry2D operator=(const StructuredGeometry2D& rhs) = delete; //Copy Assign
	virtual ~StructuredGeometry2D(); // Deconstructor

	//Output*****************************************************************
	void print() const;
	void vtkOutput(const char* fileName="output.vtk") const;
	void vtkOutput(const SpacialArray2D<double>& scalar,
		const char* fileName="output.vtk") const;
	void vtkOutput(const SpacialArray2D<double>& scalar,
		const SpacialArray2D<Vector2D<double> >& vector,
		const char* fileName="output.vtk") const;
	void vtkOutput(const SpacialArray2D<double>& scalar,
		const SpacialArray2D<Vector2D<double> >& vector,
		const SpacialArray2D<double>& divergence,
		const char* fileName="output.vtk") const;

	void vtkOutput(ofstream& output);

	//Diferential Operators**************************************************
	private:
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
	Vector2D<double> gradient(const int i, const int j,
		const SpacialArray2D<double>& scalarField) const;
	double divergence(const int i, const int j,
		const SpacialArray2D<Vector2D<double> >& vector) const;
	double divergence(const int i, const int j,
		const SpacialArray2D<double>& scalar,
		const SpacialArray2D<Vector2D<double> >& vector) const;
	Vector2D<double> divergence(const int i, const int j,
		const SpacialArray2D<double>& scalar,
		const SpacialArray2D<Vector2D<double> >& vector1,
		const SpacialArray2D<Vector2D<double> >& vector2) const;


	//Mutation Handling******************************************************
	void link_north(const StructuredGeometry2D& toLink);
	void link_south(const StructuredGeometry2D& toLink);
	void link_east(const StructuredGeometry2D& toLink);
	void link_west(const StructuredGeometry2D& toLink);

	inline void close_north()
		{ setBdy_zoneDelta_north(0.0); }

	inline void close_south()
		{ setBdy_zoneDelta_south(0.0); }

	inline void close_east()
		{ setBdy_zoneDelta_east(0.0); }

	inline void close_west()
		{ setBdy_zoneDelta_west(0.0); }

	inline void close_all_bdys()
	{
		close_north();
		close_south();
		close_east();
		close_west();
	}

	//Get Information********************************************************
	inline double get_zoneDelta_north() const
		{ return zoneDelta.dir2[numZones.get_dir2()]; }
	inline double get_zoneDelta_south() const
		{ return zoneDelta.dir2[1]; }
	inline double get_zoneDelta_east() const
		{ return zoneDelta.dir1[numZones.get_dir1()]; }
	inline double get_zoneDelta_west() const { return zoneDelta.dir1[1]; }
	inline int getCount_dir1() const { return numZones.get_dir1(); }
	inline int getCount_dir2() const { return numZones.get_dir2(); }
	inline Vector2D<int> getSize() const { return numZones; }
	inline Point2D<double> get_origin() const { return origin; }
	inline Point2D<double> get_extent() const { return extent; }

	StructuredLocalField2D<double> getLocalDeltaFactors(int i, int j) const;
	Point2D<double> getPoint(int i, int j) const;
	Vector2D<double> getMeshDelta(int i, int j) const;

	inline double get_minMeshSpacing() const { return minMeshSpacing; }

	//Coordinate System Specific Functions: Default to rectangular.
	virtual double getVolume(int i, int j) const;
	virtual StructuredLocalField2D<double> getCellAreas(int i, int j) const;

	private:
	//Set Information********************************************************
	void setBdy_zoneDelta_north(double value);
	void setBdy_zoneDelta_south(double value);
	void setBdy_zoneDelta_east(double value);
	void setBdy_zoneDelta_west(double value);

	protected:
	//Check Information******************************************************
	void checkIndices(int i, int j) const;
};

//****************************************************************Cylindrical
//****************************************************************Cylindrical
//****************************************************************Cylindrical
//StructuredCylGeometry2D*****************************************Cylindrical
/*
 * class: StructuredCylGeometry2D
 * Purpose: Hold and process geometric parameters for a structured geometry
 * 	in cylindrical coordinates.
 *
 * Description: StructuredCylGeometry2D holds information
 * 	on the spacial placement
 * 	of a particular geometry in cyilindrical. It also holds information
 * 	regarding the mesh spacing and mesh density.
 * 	Global Coordinates of structured cells
 * 	are specified using a cell-centered scheme.
 * 	This class is designed for an orthogonal grid.
 */
class StructuredCylGeometry2D : public StructuredGeometry2D
{
	public:
	StructuredCylGeometry2D(
		Point2D<double> origin=Point2D<double>(0.0,0.0),
		Point2D<double> extent=Point2D<double>(1.0,1.0),
		Vector2D<int> size=Vector2D<int>(10,10),
		Vector2D<double> refineMesh=Vector2D<double>(1.0,1.0)) 
		:	StructuredGeometry2D(origin, extent, size, refineMesh) {} // Con

	StructuredCylGeometry2D(const StructuredCylGeometry2D& that)
		:	StructuredGeometry2D(that) {}

	StructuredCylGeometry2D(const StructuredGeometry2D& that)
		:	StructuredGeometry2D(that) {}

	//Coordinate System Specific Functions.
	double getVolume(int i, int j) const;
	StructuredLocalField2D<double> getCellAreas(int i, int j) const;
};

#endif  // INCLUDE_STRUCTUREDGEOMETRY2D_H_
