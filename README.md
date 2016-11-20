# I3D08 research code

GLSL of shaders, used to implement the technique described in paper:
"Maxmimum Mipmaps for Fast, Accurate, and Scalable Dynamic Height Field
Rendering", A. Tevs, I. Ihrke, H.-P. Seidel, I3D 2008.

Also attached a simple heightmap and a colormap for testing.

Project page:
http://www.tevs.eu/project_i3d08.html


## List of files:

min_mipmap_2D_fp.glsl - computes a max-mipmap structure needed for the
algorithm. The input heightmap should have a size of (2^n+1)x(2^m+1). In
this case it is 1025x1025

mmz_vertex - is the vertex shader, which just precomputes some neccessary
values for the fragment shader

mmz_bipatch - is applied on the object in world space (NOTE: this shader
is not for the tangent space). It is applie on a box (hence 6
planes). This shader computes the intersection with the box and
initiate the algorithm

mmz_bipatch_opt_fp - is the algorithm itself incl. the intersection with
the bilinear patches.

The rest are some helping and base functions.

