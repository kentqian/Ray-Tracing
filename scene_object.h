/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		classes defining primitives in the scene

***********************************************************/

#include "util.h"
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// All primitives should provide a intersection function.
// To create more primitives, inherit from SceneObject.
// Namely, you can create, Sphere, Cylinder, etc... classes
// here.
class SceneObject {
public:
	// Returns true if an intersection occured, false otherwise.
	virtual bool intersect( Ray3D&, const Matrix4x4&, const Matrix4x4& ) = 0;
	std::vector<Point3D> Unit_AABB;
private:
	virtual void textureCoord ( Ray3D&, const Matrix4x4& ) = 0;
};

// Example primitive you can create, this is a unit square on
// the xy-plane.
class UnitSquare : public SceneObject {
public:	
	UnitSquare():t_flag(false){
		Unit_AABB.push_back(Point3D(0.5,0.5,0.0));
		Unit_AABB.push_back(Point3D(-0.5,0.5,0.0));
		Unit_AABB.push_back(Point3D(-0.5,-0.5,0.0));
		Unit_AABB.push_back(Point3D(0.5,-0.5,0.0));
	};
	UnitSquare(bool flag):t_flag(flag){
		Unit_AABB.push_back(Point3D(0.5,0.5,0.0));
		Unit_AABB.push_back(Point3D(-0.5,0.5,0.0));
		Unit_AABB.push_back(Point3D(-0.5,-0.5,0.0));
		Unit_AABB.push_back(Point3D(0.5,-0.5,0.0));
	};
	bool t_flag;
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
private:
			// Calculate texture coordinates
			void textureCoord ( Ray3D& ray, const Matrix4x4& worldToModel);
};

class UnitSphere : public SceneObject {
public:
	UnitSphere():t_flag(false){
		Unit_AABB.push_back(Point3D(-1.0,1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,-1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,-1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,-1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,-1.0,1-.0));
	};
	UnitSphere(bool flag):t_flag(flag){
		Unit_AABB.push_back(Point3D(-1.0,1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,-1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,-1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,-1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,-1.0,1-.0));
	};
	bool t_flag;
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );

private:
		// Calculate texture coordinates
		void textureCoord ( Ray3D& ray, const Matrix4x4& worldToModel);
};

class UnitCylinder : public SceneObject {
public:
	UnitCylinder():t_flag(false){
		Unit_AABB.push_back(Point3D(-1.0,1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,-1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,-1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,-1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,-1.0,1-.0));
	};
	UnitCylinder(bool flag):t_flag(flag){
		Unit_AABB.push_back(Point3D(-1.0,1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,-1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,-1.0,1.0));
		Unit_AABB.push_back(Point3D(1.0,-1.0,-1.0));
		Unit_AABB.push_back(Point3D(-1.0,-1.0,1-.0));
	};
	bool t_flag;
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );

private:
		// Calculate texture coordinates
		void textureCoord ( Ray3D& ray, const Matrix4x4& worldToModel);
};
