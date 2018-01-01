/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h,
		and the main function which specifies the
		scene to be rendered.

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <random>

bool ANTI_ALIASING = false;

Raytracer::Raytracer() : _lightSource(NULL) {
	//_root_light = new LightListNode();
	_root = new SceneDagNode();
	_tree_root = new AABBNode();
}

Raytracer::~Raytracer() {
	delete _root;
	delete _tree_root;
}

AABBNode* Raytracer::CreatingAABB(SceneDagNode* node, int flag){
    std::vector<Point3D> UnitAABB = node->obj->Unit_AABB;
	AABBNode* aabb_node = new AABBNode(node->obj,node->mat,node->trans,node->invtrans,node->modelToWorld,node->worldToModel);
	std::vector<Point3D> transformAABB;
	Matrix4x4 atrans = aabb_node->trans;
	switch(flag){
		case 0:
			// Minimimum bounding box for plane
			for (int n = 0; n < UnitAABB.size(); n++){
				transformAABB.push_back(node->trans * UnitAABB[n]);
			}
			findBounding(transformAABB,0,aabb_node);
			findBounding(transformAABB,1,aabb_node);
			findBounding(transformAABB,2,aabb_node);
			calculatesurfacearea(aabb_node);
			return aabb_node;
		break;
		case 1:
			// Minimimum bounding box for sphere
			aabb_node->minX = atrans[0][3] - sqrt(pow(atrans[0][0],2.0) + pow(atrans[0][1],2.0) + pow(atrans[0][2],2.0));
			aabb_node->maxX = atrans[0][3] + sqrt(pow(atrans[0][0],2.0) + pow(atrans[0][1],2.0) + pow(atrans[0][2],2.0));
			aabb_node->minY = atrans[1][3] - sqrt(pow(atrans[1][0],2.0) + pow(atrans[1][1],2.0) + pow(atrans[1][2],2.0));
			aabb_node->maxY = atrans[1][3] + sqrt(pow(atrans[1][0],2.0) + pow(atrans[1][1],2.0) + pow(atrans[1][2],2.0));
			aabb_node->minZ = atrans[2][3] - sqrt(pow(atrans[2][0],2.0) + pow(atrans[2][1],2.0) + pow(atrans[2][2],2.0));
			aabb_node->maxZ = atrans[2][3] + sqrt(pow(atrans[2][0],2.0) + pow(atrans[2][1],2.0) + pow(atrans[2][2],2.0));
			calculatesurfacearea(aabb_node);
			return aabb_node;
		break;
	}
}

void Raytracer::findBounding(std::vector<Point3D> transformAABB, int pos, AABBNode* aabb_node){
	// Fill the minimimum bounding values when calculating plane bounding box
	double temp_min = transformAABB.front()[pos];
	double temp_max = transformAABB.front()[pos];

	for(int n = 0; n < transformAABB.size(); n++){
		if(transformAABB[n][pos] < temp_min){
			temp_min = transformAABB[n][pos];
		}
		if(transformAABB[n][pos] > temp_max){
			temp_max = transformAABB[n][pos];
		}
	}
	switch(pos){
		case 0:
			aabb_node->minX = temp_min;
			aabb_node->maxX = temp_max;
		break;
		case 1:
			aabb_node->minY = temp_min;
			aabb_node->maxY = temp_max;
		break;
		case 2:
			aabb_node->minZ = temp_min;
			aabb_node->maxZ = temp_max;
		break;
	}
}

void Raytracer::AABBinsert(AABBNode* leafNode){
	if(!has_root){
		// for tree root
		has_root = true;
		_tree_root = leafNode;
		return;
	}
	AABBNode* treeRoot = _tree_root;

	// Search for the best place to put the new leaf in the tree
	// We use surface area and depth as search heuristics
	while(!treeRoot->isLeaf()){
		// Because of the test in the while loop above we know we are never a leaf inside it
		AABBNode* lefNode = treeRoot->lef_child;
		AABBNode* rigNode = treeRoot->rig_child;

		AABBNode* combinedAabb = AABBMerge(treeRoot, leafNode);

		double newParentNodeCost = 2.0 * combinedAabb->surface_area;
		double minimumPushDownCost = 2.0 * (combinedAabb->surface_area - treeRoot->surface_area);
		delete combinedAabb;

		// Use the costs to figure out whether to create a new parent here or descend
		double costLeft;
		double costRight;
		AABBNode* lefMerge = AABBMerge(leafNode, lefNode);
		AABBNode* rigMerge = AABBMerge(leafNode, rigNode);

		if(lefNode->isLeaf()){
			costLeft = lefMerge->surface_area + minimumPushDownCost;
		}else{
			costLeft = (lefMerge->surface_area - lefNode->surface_area) + minimumPushDownCost;
		}
		if(rigNode->isLeaf()){
			costRight = rigMerge->surface_area + minimumPushDownCost;
		}else{
			costRight = (lefMerge->surface_area - lefNode->surface_area) + minimumPushDownCost;
		}

		// If the cost of creating a new parent node here is less than descending in either direction then
		// we know we need to create a new parent node, errrr, here and attach the leaf to that
		if(newParentNodeCost < costLeft && newParentNodeCost < costRight){
			break;
		}
		delete lefMerge;
		delete rigMerge;

		// Otherwise descend in the cheapest direction
		if(costLeft < costRight){
			treeRoot = lefNode;
		}else{
			treeRoot = rigNode;
		}
	}

	// The leafs sibling is going to be the node we found above and we are going to create a new
	// parent node and attach the leaf and this item
	AABBNode* leafSibling = treeRoot;
	AABBNode* oldParent = leafSibling->parent;
	AABBNode* newParent = AABBMerge(leafSibling,leafNode);
	newParent->parent = oldParent;
	newParent->lef_child = leafSibling;
	newParent->rig_child = leafNode;
	leafNode->parent = newParent;
	leafSibling->parent = newParent;

	if(oldParent == NULL){

		// the old parent was the root and so this is now the root
		_tree_root = newParent;
	}else{

		// The old parent was not the root and so we need to patch the left or right index to
		// point to the new node
		if(oldParent->lef_child == leafSibling){
			oldParent->lef_child = newParent;
		}else{
			oldParent->rig_child = newParent;
		}
	}
	// Finally we need to walk back up the tree fixing heights and areas
	fixUpwardsTree(newParent);
}

void Raytracer::fixUpwardsTree(AABBNode* treeNode){
	while(treeNode != NULL){
		// Fix height and area
		AABBNode* lefNode = treeNode->lef_child;
		AABBNode* rigNode = treeNode->rig_child;
		AABBNode* updateNode = AABBMerge(lefNode,rigNode);
		treeNode->minX = updateNode->minX;
		treeNode->minY = updateNode->minY;
		treeNode->minZ = updateNode->minZ;
		treeNode->maxX = updateNode->maxX;
		treeNode->maxY = updateNode->maxY;
		treeNode->maxZ = updateNode->maxZ;
		treeNode->surface_area = updateNode->surface_area;
		treeNode = treeNode->parent;
	}
}

AABBNode* Raytracer::AABBMerge(AABBNode* first, AABBNode* second){
	double minX = std::min(first->minX, second->minX);
	double minY = std::min(first->minY, second->minY);
	double minZ = std::min(first->minZ, second->minZ);
	double maxX = std::max(first->maxX, second->maxX);
	double maxY = std::max(first->maxY, second->maxY);
	double maxZ = std::max(first->maxZ, second->maxZ);
	AABBNode* newNode = new AABBNode(minX,minY,minZ,maxX,maxY,maxZ);
	calculatesurfacearea(newNode);
	return newNode;
}

void Raytracer::calculatesurfacearea(AABBNode* aabb_node){
	double width = aabb_node->maxX - aabb_node->minX;
	double height = aabb_node->maxY - aabb_node->minY;
	double depth = aabb_node->maxZ - aabb_node->minZ;
	aabb_node->surface_area =  2.0 * (width * height + width * depth + height * depth);
}

bool Raytracer::intersectBB(AABBNode* aabb, Ray3D& ray) {
	// Reference:
	//https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection

	Point3D tmin, tmax;
	Point3D invdir(1/ray.dir[0], 1/ray.dir[1], 1/ray.dir[2]);
	Point3D aabbMin(aabb->minX, aabb->minY, aabb->minZ);
	Point3D aabbMax(aabb->maxX, aabb->maxY, aabb->maxZ);

	for (int i = 0; i < 3; i++){
		if (invdir[i] >= 0) {
			tmin[i] = (aabbMin[i] - ray.origin[i]) * invdir[i];
			tmax[i] = (aabbMax[i] - ray.origin[i]) * invdir[i];
		} else {
			tmin[i] = (aabbMax[i] - ray.origin[i]) * invdir[i];
			tmax[i] = (aabbMin[i] - ray.origin[i]) * invdir[i];
		}
	}

	if ((tmin[0] > tmax[1]) || (tmin[1] > tmax[0])) {
		return false;
	}
	if (tmin[1] > tmin[0]) {
		tmin[0] = tmin[1];
	}
	if (tmax[1] < tmax[0]) {
		tmax[0] = tmax[1];
	}

	if ((tmin[0] > tmax[2]) || (tmin[2] > tmax[0])) {
		return false;
	}
	if (tmin[2] > tmin[0]) {
		tmin[0] = tmin[2];
	}
	if (tmax[2] < tmax[0]) {
		tmax[0] = tmax[2];
	}

	return true;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent,
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;

	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}

	return node;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;

	if ( _lightSource == NULL) {
		_lightSource = new LightListNode( light, tmp );
	}else {
		while (_lightSource->next != NULL) {
			_lightSource = _lightSource->next;
		}
		_lightSource = new LightListNode( light, tmp);
	}
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;

	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation;
			angle = -angle;
		}
		else {
			node->invtrans = rotation*node->invtrans;
		}
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;

	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation;
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans;
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;

	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale;
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans;
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view,
		Vector3D up ) {
	Matrix4x4 mat;
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat;
}


void Raytracer::computeTransforms( SceneDagNode* node )
{
    SceneDagNode *childPtr;
    if (node->parent != NULL )
    {
        node->modelToWorld = node->parent->modelToWorld*node->trans;
        node->worldToModel = node->invtrans*node->parent->worldToModel;
    }
    else
    {
        node->modelToWorld = node->trans;
        node->worldToModel = node->invtrans;
    }
    // Traverse the children.
    childPtr = node->child;
    while (childPtr != NULL) {
        computeTransforms(childPtr);
        childPtr = childPtr->next;
    }

}

void Raytracer::computeTransformsBB( AABBNode* node ){
	if (node->parent != NULL){
		node->modelToWorld = node->parent->modelToWorld*node->trans;
        node->worldToModel = node->invtrans*node->parent->worldToModel;
	}
	else{
		node->modelToWorld = node->trans;
        node->worldToModel = node->invtrans;
	}
	AABBNode* lef_node = node->lef_child;
	AABBNode* rig_node = node->rig_child;
	if (lef_node != NULL){
		computeTransformsBB(lef_node);
		computeTransformsBB(rig_node);
	}
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
    SceneDagNode *childPtr;

    // Applies transformation of the current node to the global
    // transformation matrices.
    if (node->obj) {
        // Perform intersection.
        if (node->obj->intersect(ray, node->worldToModel, node->modelToWorld)) {
            ray.intersection.mat = node->mat;
        }
    }
    // Traverse the children.
    childPtr = node->child;
    while (childPtr != NULL) {
        traverseScene(childPtr, ray);
        childPtr = childPtr->next;
    }

}

void Raytracer::traverseSceneBB( AABBNode* aabb, Ray3D& ray ) {
	if (aabb != NULL) {
		if(intersectBB(aabb, ray)) {
			// if intersect with aabb, find recursively
			traverseSceneBB(aabb->lef_child, ray);
			traverseSceneBB(aabb->rig_child, ray);

			if (aabb->obj) {
				if (aabb->obj->intersect(ray, aabb->worldToModel, aabb->modelToWorld)) {
					ray.intersection.mat = aabb->mat;
				}
			}
		}
	}
}

void Raytracer::computeShading( Ray3D& ray ) {
    LightListNode* curLight = _lightSource;

    for (;;) {
        if (curLight == NULL) break;
        // Each lightSource provides its own shading function.
        // Implement shadows here if needed.

		// Implement soft shadow, sampling the light direction with ramdom jitter
		Colour colour(0.0, 0.0, 0.0);
		for (float i = -1; i < 1.0; i = i + 0.1) {
			// double x = rand() / (float) RAND_MAX - 0.5;
			// double y = rand() / (float) RAND_MAX - 0.5;
			// double z = rand() / (float) RAND_MAX - 0.5;
			// Vector3D random(x, y, z);

			Vector3D shadow_dir = curLight->light->get_position() - ray.intersection.point;
			for (int j=0; j<3; j++) {shadow_dir[j] += i;}
			shadow_dir.normalize();

			Point3D shadow_point = ray.intersection.point + 0.0001 * shadow_dir;
			Ray3D shadowRay(shadow_point,shadow_dir);

			// traverseScene(_root,shadowRay);

			// AABB traversal
			traverseSceneBB(_tree_root,shadowRay);

			curLight->light->shade(ray);

			if (shadowRay.intersection.none) {
				colour = colour + (0.05 * ray.col);
			}
		}
		ray.col = colour;

		// Without shadow feature(Only Phong and reflection)
		// curLight->light->shade(ray);

        curLight = curLight->next;
    }
}

void Raytracer::initPixelBuffer() {
    int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
    _rbuffer = new unsigned char[numbytes];
    std::fill_n(_rbuffer, numbytes,0);
    _gbuffer = new unsigned char[numbytes];
    std::fill_n(_gbuffer, numbytes,0);
    _bbuffer = new unsigned char[numbytes];
    std::fill_n(_bbuffer, numbytes,0);
}

void Raytracer::flushPixelBuffer( std::string file_name ) {
    bmp_write( file_name.c_str(), _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
    delete _rbuffer;
    delete _gbuffer;
    delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray, int recur_time) {
    Colour col(0.0, 0.0, 0.0);
    Colour reflection_col(0.0, 0.0, 0.0);

    // traverseScene(_root, ray);

	// AABB traversal
	traverseSceneBB(_tree_root,ray);

	// recur_time controls the reflection times
    if (ray.intersection.none || recur_time >= 3){
    	return col;
    }

		// Add texture color if texture flag is true
		Colour texture_col = Colour (1, 1, 1);
		if (ray.intersection.texture_flag) {
			texture_col = get_texture_color(ray);
		}

    // Don't bother shading if the ray didn't hit
    // anything.

    if (!ray.intersection.none) {
        computeShading(ray);

        Vector3D N = ray.intersection.normal;
        Vector3D D = ray.dir;
        Vector3D R = D - (2 * N.dot(D)) * N;
        N.normalize();
        D.normalize();
        R.normalize();

        Point3D R_org = ray.intersection.point + 0.0001 * R;

        Ray3D R_ray(R_org,R);

        if (N.dot(R) > 0){
        	reflection_col = pow(ray.intersection.mat->reflect_ratio, recur_time + 1) * shadeRay(R_ray, recur_time + 1);
    	}
    }

	// Add texture color and turn down the brightness
	if (ray.intersection.texture_flag) {
		col = 0.6 * ray.col + reflection_col + 0.8 * texture_col;
	} else {
		col = ray.col + reflection_col;
	}


    col.clamp();
    return col;
}

void Raytracer::read_texture(char texture_name[]) {
		// From TA: https://bb-2017-09.teach.cs.toronto.edu/t/loading-textures-models-for-a3/1891
		rarray = NULL;
		garray = NULL;
		barray = NULL;

		if (bmp_read(texture_name, &t_width,&t_height,&rarray,&garray,&barray)) {
			printf("Error reading the texture image");
		}
}

Colour Raytracer::get_texture_color(Ray3D& ray) {

		if (ray.intersection.uvCoord[0] && ray.intersection.uvCoord[1]) {
			// Get x, y coordinate of pixel
		  int x = (ray.intersection.uvCoord[0]) * t_width;
			int y = (ray.intersection.uvCoord[1]) * t_height;

			// Index of rgb arrays
			int i = x * t_width + y;

			// Get r, g, b texture color at index i
			double r = (1.0/255) * rarray[i];
		  double g = (1.0/255) * garray[i];
			double b = (1.0/255) * barray[i];

			return Colour(r, g, b);
		} else {
			return Colour(0, 0, 0);
		}

}

void Raytracer::render( int width, int height, Point3D eye, Vector3D view,
        Vector3D up, double fov, std::string fileName ) {
    // computeTransforms(_root);
	computeTransformsBB(_tree_root);
    Matrix4x4 viewToWorld;
    _scrWidth = width;
    _scrHeight = height;
    double factor = (double(height)/2)/tan(fov*M_PI/360.0);

		read_texture("blue.bmp");

    initPixelBuffer();
    viewToWorld = initInvViewMatrix(eye, view, up);

    // Construct a ray for each pixel.
    for (int i = 0; i < _scrHeight; i++) {
        for (int j = 0; j < _scrWidth; j++) {
            // Sets up ray origin and direction in view space,
            // image plane is at z = -1.
						Point3D origin(0, 0, 0);
						Point3D imagePlane;
						Colour col;
						//int anti = ANTI_ALIASING;
						//double co_anti;

						// Implementation of anti-aliasing by using the super sampling method, references： https://www.ics.uci.edu/~gopi/CS211B/RayTracing%20tutorial.pdf Page 12
						if (ANTI_ALIASING) {
							for (float m = i; m < i + 1.0f; m += 0.5f) {
								for (float n = j; n < j + 1.0f; n += 0.5f) {
									//co_anti = 1/pow(anti, 2);

									imagePlane[0] = (-double(width)/2 + 0.5 + n)/factor;
									imagePlane[1] = (-double(height)/2 + 0.5 + m )/factor;
									imagePlane[2] = -1;

									// TODO: Convert ray to world space and call
									// shadeRay(ray) to generate pixel colour.

									Ray3D ray;

									ray.dir =  viewToWorld * (imagePlane - origin);
									ray.origin = viewToWorld * origin;
									ray.dir.normalize();

									col = shadeRay(ray, 0);

									_rbuffer[i*width+j] += int(col[0]*255*0.25f);
									_gbuffer[i*width+j] += int(col[1]*255*0.25f);
									_bbuffer[i*width+j] += int(col[2]*255*0.25f);

								}
							}
						} else {
							imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
							imagePlane[1] = (-double(height)/2 + 0.5 + i )/factor;
							imagePlane[2] = -1;

							// TODO: Convert ray to world space and call
							// shadeRay(ray) to generate pixel colour.

							Ray3D ray;

							ray.dir =  viewToWorld * (imagePlane - origin);
							ray.origin = viewToWorld * origin;
							ray.dir.normalize();

							col = shadeRay(ray, 0);
							_rbuffer[i*width+j] = int(col[0]*255);
							_gbuffer[i*width+j] = int(col[1]*255);
							_bbuffer[i*width+j] = int(col[2]*255);
						}
		}
	}

	flushPixelBuffer(fileName);
}
