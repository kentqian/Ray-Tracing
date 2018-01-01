/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0),
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point,
	// intersection.normal, intersection.none, intersection.t_value.
	//
	// HINT: Remember to first transform the ray into object space
	// to simplify the intersection test.

	// Define ray's origin and direction in the object space

	Ray3D rayObject;
	rayObject.origin = worldToModel * ray.origin;
	rayObject.dir = worldToModel * ray.dir;

	double t = - rayObject.origin[2] / rayObject.dir[2];

	// If ray intersect with the UnitSquare
	if (t > 0) {
		Point3D p;
		p = rayObject.origin + t * rayObject.dir;

		if (p[0] >= -0.5 && p[0] <= 0.5 && p[1] >= -0.5 && p[1] <= 0.5) {
			if(ray.intersection.none || t < ray.intersection.t_value) {

				// Implement intersection code for UnitSquare
				ray.intersection.point = modelToWorld * p;
				ray.intersection.t_value = t;
				ray.intersection.none = false;

				Vector3D n(0.0, 0.0, 1.0);
				ray.intersection.normal = worldToModel.transpose() * n;
				ray.intersection.normal.normalize();
				if (t_flag) {
					textureCoord(ray, worldToModel);
					ray.intersection.texture_flag = true;
				}
				return true;
			}
		}
	}
	return false;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred
	// on the origin.
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point,
	// intersection.normal, intersection.none, intersection.t_value.
	//
	// HINT: Remember to first transform the ray into object space
	// to simplify the intersection test.

	Ray3D rayObject;
	rayObject.origin = worldToModel * ray.origin;
	rayObject.dir = worldToModel * ray.dir;

	Point3D sphereOrigin(0, 0, 0);
	Vector3D distance = rayObject.origin - sphereOrigin;

	double a;
	double b;
	double c;
	double delta;
	double t;

	// Determine how many ts depend on delta
	a = rayObject.dir.dot(rayObject.dir);
	b = 2 * rayObject.dir.dot(distance);
	c = distance.dot(distance) - 1;
	delta = pow(b,2) - 4*a*c;

	// If ray intersect with the UnitSphere
	if(delta == 0){

		// Only one t
		t = -b/(2*a);
	}else if(delta > 0){

		// Select minimum valid t from two ts
		double t0 = (-b - sqrt(delta))/(2*a);
		double t1 = (-b + sqrt(delta))/(2*a);
		if(t0 > 0){
			t = t0;
		}else if(t0 < 0 && t1 > 0){
			t = t1;
		}else{
			return false;
		}

	}else{
		return false;
	}

	if(ray.intersection.none || t < ray.intersection.t_value){

		// Implement intersection code for UnitSphere
		ray.intersection.none = false;
		ray.intersection.t_value = t;
		Point3D p = rayObject.origin+ray.intersection.t_value * rayObject.dir;
		ray.intersection.point = modelToWorld * p;
		ray.intersection.normal = worldToModel.transpose() * (p - sphereOrigin);
		ray.intersection.normal.normalize();
		if (t_flag) {
			textureCoord(ray, worldToModel);
			ray.intersection.texture_flag = true;
		}
		return true;
	}
	return false;

}

bool UnitCylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel, const Matrix4x4& modelToWorld ){
	Ray3D rayObject;
	rayObject.origin = worldToModel * ray.origin;
	rayObject.dir = worldToModel * ray.dir;
	
	// Firstly, Find t0 and t1 whether the intersection crosses with an infinite cylinder
	double a, b, c, delta, t0, t1;

	a = rayObject.dir[0]*rayObject.dir[0] + rayObject.dir[2]*rayObject.dir[2];
	b = 2 * rayObject.origin[0]*rayObject.dir[0] + 2 * rayObject.origin[2] * rayObject.dir[2];
	c = rayObject.origin[0]*rayObject.origin[0] + rayObject.origin[2]*rayObject.origin[2] - 1;

	delta = pow(b,2.0) - 4*a*c;
	if(delta < 0){
		return false;
		}

	t0 = (-b + sqrt(delta)) / (2*a);
	t1 = (-b - sqrt(delta)) / (2*a);
	if(t0 > t1){
		double tmp = t0;
		t0=t1;
		t1=tmp;
	}

	// Secondly, Find y0 and y1 as the cylinder cap situation
	double y0 = rayObject.origin[1] + t0*rayObject.dir[1];
	double y1 = rayObject.origin[1] + t1*rayObject.dir[1];

	if(y0 < -1){
		if(y1 < -1){
			
			// Miss down cap 
			return false;
		}else{
			double t = t0 + (t1-t0) * (y0+1) / (y0-y1);
			if(t <= 0){
				return false;
			}
			if(ray.intersection.none || t < ray.intersection.t_value){
				ray.intersection.none = false;
				ray.intersection.t_value = t;
				Point3D p = rayObject.origin + ray.intersection.t_value * rayObject.dir;
				ray.intersection.point = modelToWorld * p;
				ray.intersection.normal = worldToModel.transpose() * Vector3D(0.0, -1.0, 0.0);
				ray.intersection.normal.normalize();
				if (t_flag) {
					textureCoord(ray, worldToModel);
					ray.intersection.texture_flag = true;
				}
				return true;
			}
		}
	}else if(y0 > -1 && y0 < 1){
		
		// Hit the side of cylinder
		if(t0 <= 0){
			return false;
		}
		if(ray.intersection.none || t0 < ray.intersection.t_value){
			ray.intersection.none = false;
			ray.intersection.t_value = t0;
			Point3D p = rayObject.origin + ray.intersection.t_value * rayObject.dir;
			ray.intersection.point = modelToWorld * p;
			ray.intersection.normal = worldToModel.transpose() * Vector3D(p[0],0.0,p[2]);
			ray.intersection.normal.normalize();
			if (t_flag) {
				textureCoord(ray, worldToModel);
				ray.intersection.texture_flag = true;
			}
			return true;
		}
	}
	else if(y0 > 1){
		if(y1 > 1){
			// Miss top cap
			return false;
		}else{
			double t = t0 + (t1-t0) * (y0-1) / (y0-y1);
			if(t <= 0){
				return false;
			}
			if(ray.intersection.none || t < ray.intersection.t_value){
				ray.intersection.none = false;
				ray.intersection.t_value = t;
				Point3D p = rayObject.origin + ray.intersection.t_value * rayObject.dir;
				ray.intersection.point = modelToWorld * p;
				ray.intersection.normal = worldToModel.transpose() * Vector3D(0.0, 1.0, 0.0);
				ray.intersection.normal.normalize();
				if (t_flag) {
					textureCoord(ray, worldToModel);
					ray.intersection.texture_flag = true;
				}
				return true;
			}
		}
	}
	return false;
}

void UnitSphere::textureCoord(Ray3D& ray, const Matrix4x4& worldToModel) {
	// Calculation reference: https://www.cs.unc.edu/~rademach/xroads-RT/RTarticle.html#texturemap
	Vector3D normal = ray.intersection.normal;
	Vector3D v_n(0, 1, 0);
	Vector3D v_e(0, 0, 1);
	Vector3D v_p(normal);

	double u, v;
	double phi = acos(-((v_n).dot(v_p)));
	v = phi / M_PI;

	double theta = (acos(((v_p).dot(v_e))/sin(phi))) / (2 * M_PI);
	if (v_p.dot( v_n.cross(v_e)) > 0) {
		u = theta;
	} else {
		u = 1 - theta;
	}
	ray.intersection.uvCoord = Point3D(u, v, 0);
}

void UnitSquare::textureCoord(Ray3D& ray, const Matrix4x4& worldToModel) {
	Point3D point = worldToModel * ray.intersection.point;
	double u = point[0] + 0.5;
	double v = 0.5 - point[1];
	ray.intersection.uvCoord = Point3D(u, v, 0);
}

void UnitCylinder::textureCoord(Ray3D& ray, const Matrix4x4& worldToModel) {
	Vector3D normal = ray.intersection.normal;
	Point3D point = worldToModel * ray.intersection.point;
	Vector3D v_e(0, 0, 1);

	double v, theta;
	theta = acos(v_e.dot(normal)) / (2 * M_PI);
	v = point[1];
	ray.intersection.uvCoord = Point3D(theta, v, 0);
}
