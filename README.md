# Introduction:

This is an advanced [ray tracing] demo based on Phong's shading, implementing various physical features, such as reflection, compound object(cylinders, balls, planes), anti-aliasing, soft shadows and texture mapping.

#  File Structure:

(1) Cylinder: Done in scene_object.cpp, Unitcylinder() function

(2) Anti-Aliasing: Done in raytracer.cpp, render() function  

(3) Soft Shadow: Done in raytracer.cpp, computeShading() function
 
(4) Texture Mapping: 
	- Texture Coordinates are calculated in scene_object.cpp
	- Texture Colour is calculated in get_texture_color in raytracer.cpp
	- Read texture file is done with a helper function read_texture() in raytracer.cpp 
	- Shade Texture color is done in function shade_ray() in raytracer.cpp 

(5) BSP tree: Done in raytracer.cpp, traverseSceneBB() function

# File Descriptions:

raytracer.{cpp, h} 
The main body of the raytracer, including the scene graph. 

util.{cpp, h}
Includes definitions and implementations for points, vectors, matrices, 
and some structures for intersection and shading.  

light_source.{cpp, h}
Defines the basic light class. It defines different types of 
lights, which shades the ray differently.  Point lights are sufficient 
for most scenes.  

scene_object.{cpp, h}
Defines object primitives in the scene 
Implements the intersect function which finds the intersection point 
between the ray and the primitive. 

bmp_io.{cpp, h}
I/O functions for .bmp files.

# References:

(1) BSP tree: (Traversal function does not work very well, while AABB tree was completely implemented)
- AABB tree is based on the surface area and depth of node as search heuristics. 
- Reference: https://github.com/JamesRandall/SimpleVoxelEngine/blob/master/voxelEngine/src/AABBTree.cpp

(2) Cylinder: 
- Reference: http://woo4.me/wootracer/cylinder-intersection/

(3) Anti-Aliasing:  
- Anti-Aliasing is implemented by using the supersampling method. 
- To see the Anti-Aliasing effect, change the boolean ANTI_ALIASING to true (in raytracer.cpp on the top). 
- Reference: https://www.ics.uci.edu/~gopi/CS211B/RayTracing%20tutorial.pdf (Page 12)

(4) Soft Shadow
- Soft Shadow is implemented by sampling the light direction with random jitter

(5) Texture Mapping
- A square and sphere texture mappings are implemented, and each is controlled by a flag in Intersection struct, hence, you can see both texture mapping by setting until.n file, texture_flag to true, and square_texture or sphere_texture to true depends on what you want to map. Setting texture_flag to false will disable the texture mapping for all objects. 
- The overall method to produce texture mapping is to map every point on an object to a pixel on the image, then get the colour of the pixel and shade. 
	- Reference: https://www.cs.unc.edu/~rademach/xroads-RT/RTarticle.html#texturemap


[ray tracing]: <https://en.wikipedia.org/wiki/Ray_tracing_(graphics)>