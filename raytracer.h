/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		This file contains the interface and
		datastructures of the raytracer.
		Simple traversal and addition code to
		the datastructures are given to you.

***********************************************************/

#include "util.h"
#include "scene_object.h"
#include "light_source.h"


// Linked list containing light sources in the scene.
struct LightListNode {
	LightListNode() : light(NULL), next(NULL) {}
	LightListNode( LightSource* light, LightListNode* next = NULL ) :
		light(light), next(next) {}
	~LightListNode() {
		if (!light) delete light;
	}
	LightSource* light;
	LightListNode* next;
};

// The scene graph, containing objects in the scene.
struct SceneDagNode {
	SceneDagNode() :
		obj(NULL), mat(NULL),
		next(NULL), parent(NULL), child(NULL) {
	}

	SceneDagNode( SceneObject* obj, Material* mat ) :
		obj(obj), mat(mat), next(NULL), parent(NULL), child(NULL) {
		}

	~SceneDagNode() {
		if (!obj) delete obj;
		if (!mat) delete mat;
	}

	// Pointer to geometry primitive, used for intersection.
	SceneObject* obj;
	// Pointer to material of the object, used in shading.
	Material* mat;
	// Each node maintains a transformation matrix, which maps the
	// geometry from object space to world space and the inverse.
	Matrix4x4 trans;
	Matrix4x4 invtrans;
	Matrix4x4 modelToWorld;
	Matrix4x4 worldToModel;

	// Internal structure of the tree, you shouldn't have to worry
	// about them.
	SceneDagNode* next;
	SceneDagNode* parent;
	SceneDagNode* child;
};

// The AABB Node
struct AABBNode {
public:
	// Pointer to geometry primitive, used for intersection.
	SceneObject* obj;
	// Pointer to material of the object, used in shading.
	Material* mat;
	// Each node maintains a transformation matrix, which maps the 
	// geometry from object space to world space and the inverse.
	Matrix4x4 trans;
	Matrix4x4 invtrans;
	Matrix4x4 modelToWorld;
	Matrix4x4 worldToModel;
	// Bounding value
	double minX;
	double minY;
	double minZ;
	double maxX;
	double maxY;
	double maxZ;
	double surface_area;

	AABBNode* parent;
	AABBNode* lef_child;
	AABBNode* rig_child; 

	AABBNode() 
		: obj(NULL), mat(NULL), trans(), invtrans(), modelToWorld(), worldToModel(),
		minX(0.0), minY(0.0), minZ(0.0), maxX(0.0), maxY(0.0), maxZ(0.0),surface_area(0.0),
		parent(NULL),lef_child(NULL),rig_child(NULL){};

	AABBNode( SceneObject* obj, Material* mat, Matrix4x4 trans, Matrix4x4 invtrans, Matrix4x4 modelToWorld, Matrix4x4 worldToModel) : 
	obj(obj), mat(mat), trans(trans), invtrans(invtrans), modelToWorld(modelToWorld), worldToModel(worldToModel),
	minX(0.0), minY(0.0), minZ(0.0), maxX(0.0), maxY(0.0), maxZ(0.0),surface_area(0.0),
	parent(NULL),lef_child(NULL),rig_child(NULL){};
	
	AABBNode(double minX, double minY, double minZ, double maxX, double maxY, double maxZ) :
	obj(NULL), mat(NULL), trans(), invtrans(), modelToWorld(), worldToModel(),
	minX(minX), minY(minY), minZ(minZ), maxX(maxX), maxY(maxY), maxZ(maxZ),surface_area(0.0),
	parent(NULL),lef_child(NULL),rig_child(NULL){};

	bool isLeaf() const { return lef_child == NULL; }

	bool operator==(AABBNode* other) const {
		return (obj == other->obj);
	}
};

class Raytracer {
public:
	Raytracer();
	~Raytracer();

	// Renders an image fileName with width and height and a camera
	// positioned at eye, with view vector view, up vector up, and
	// field of view fov.
	void render( int width, int height, Point3D eye, Vector3D view,
			Vector3D up, double fov, std::string fileName );

	// Add an object into the scene, with material mat.  The function
	// returns a handle to the object node you just added, use the
	// handle to apply transformations to the object.
	SceneDagNode* addObject( SceneObject* obj, Material* mat ) {
		return addObject(_root, obj, mat);
	}

	// Add an object into the scene with a specific parent node,
	// don't worry about this unless you want to do hierarchical
	// modeling.  You could create nodes with NULL obj and mat,
	// in which case they just represent transformations.
	SceneDagNode* addObject( SceneDagNode* parent, SceneObject* obj,
			Material* mat );

	// LightListNode* addLightSource( LightSource* light) {
	// 	return addLightSource(_root_light, light);
	// }

	// Add a light source.
	LightListNode* addLightSource( LightSource* light );

	// Transformation functions are implemented by right-multiplying
	// the transformation matrix to the node's transformation matrix.

	// Apply rotation about axis 'x', 'y', 'z' angle degrees to node.
	void rotate( SceneDagNode* node, char axis, double angle );

	// Apply translation in the direction of trans to node.
	void translate( SceneDagNode* node, Vector3D trans );

	// Apply scaling about a fixed point origin.
	void scale( SceneDagNode* node, Point3D origin, double factor[3] );

	// Read the texture file with name texture_name
	void read_texture(char texture_name[]);

	// Find the bounding values of transformed Object
	void findBounding(std::vector<Point3D> transformAABB, int pos, AABBNode* aabb_node);

	// Create a new AABBnode
	AABBNode* CreatingAABB(SceneDagNode* node, int flag);

	// Calculate the surface area of AABBnode
	void calculatesurfacearea(AABBNode* aabb_node);

	// Check whether hit the Bounding box
	bool intersectBB(AABBNode * aabb, Ray3D & ray);

	// Merge two AABBnode to a new AABBnode
	AABBNode* AABBMerge(AABBNode* first, AABBNode* second);

	// Insert AABBnode to AABBtree
	void AABBinsert(AABBNode* leafNode);

	// Update the surface area and related data upwards to the AABBtree root
	void fixUpwardsTree(AABBNode* treeNode);

	// traverseScene with AABB
	void traverseSceneBB(AABBNode * aabb, Ray3D & ray);

	void computeTransformsBB( AABBNode* node );
	
	AABBNode *_tree_root;


private:
	// Allocates and initializes the pixel buffer for rendering, you
	// could add an interesting background to your scene by modifying
	// this function.
	void initPixelBuffer();

	// Saves the pixel buffer to a file and deletes the buffer.
	void flushPixelBuffer(std::string file_name);

	// Return the colour of the ray after intersection and shading, call
	// this function recursively for reflection and refraction.
	Colour shadeRay( Ray3D& ray, int recur_time);

	// Constructs a view to world transformation matrix based on the
	// camera parameters.
	Matrix4x4 initInvViewMatrix( Point3D eye, Vector3D view, Vector3D up );

	// Traversal code for the scene graph, the ray is transformed into
	// the object space of each node where intersection is performed.
	void traverseScene( SceneDagNode* node, Ray3D& ray );

	// After intersection, calculate the colour of the ray by shading it
	// with all light sources in the scene.
    void computeShading( Ray3D& ray );

    // Precompute the modelToWorld and worldToModel transformations for each
    // object in the scene.
    void computeTransforms( SceneDagNode* node );


    // Width and height of the viewport.
    int _scrWidth;
    int _scrHeight;

    // Light list and scene graph.
    LightListNode *_lightSource;
    SceneDagNode *_root;

    // Pixel buffer.
    unsigned char* _rbuffer;
    unsigned char* _gbuffer;
    unsigned char* _bbuffer;

	// Texture pixel rgb arrays
	unsigned char* rarray;
	unsigned char* garray;
	unsigned char* barray;

	// Texture map width and height
	unsigned long int t_width;
	long int t_height;

    // Maintain global transformation matrices similar to OpenGL's matrix
    // stack.  These are used during scene traversal.
    Matrix4x4 _modelToWorld;
    Matrix4x4 _worldToModel;

	// Anti-aliasing number
	int anti;

	//AABB root flag
	bool has_root = false;

	// Get the texture color at ray coordinates
	Colour get_texture_color(Ray3D& ray);
};
