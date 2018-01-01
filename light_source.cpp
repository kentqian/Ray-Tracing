/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"
#include <algorithm>

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray
	// is available.  So be sure that traverseScene() is called on the ray
	// before this function.

	Vector3D N = ray.intersection.normal;
	Vector3D V = -ray.dir;
	Vector3D L = get_position() - ray.intersection.point;
	Vector3D H = L + V;

	//normalize each vector function need
	N.normalize();
	V.normalize();
	L.normalize();
	H.normalize();

	// each part of phong shading implementation
	Colour ambient = ray.intersection.mat->ambient * _col_ambient;
	Colour diffuse = ray.intersection.mat->diffuse * (std::max(0.0,N.dot(L)) * _col_diffuse);
	Colour specular = ray.intersection.mat->specular * (std::max(0.0,pow(N.dot(H),ray.intersection.mat->specular_exp)) * _col_specular);

	// render with scene signuature
	// ray.col = ray.intersection.mat->diffuse;
	
	// render colour without specular
	// ray.col = ambient + diffuse;

	// render colour with phong shading
	ray.col = ambient + diffuse + specular;
	
	ray.col.clamp();

}
