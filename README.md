# Computer Graphics: Virtual Camera
## Intro
This project is created for the university class named Computer Graphics at Aristotle University of Thessaloniki (AUTh). It's the second out of three repositories referenced on the same class.

## General
The PYTHON folder that is uploaded contains 3 Python files:
* Demo.py
* Camera_functions.py
* Triangle_filling.py
The first file, Demo.py, is the one that calls the functions to perform the requested operations. Initially, it reads the data from the available .npy file and stores it in the corresponding variables. It also initializes certain constants defined in the assignment. Immediately afterward, it captures images of the object successively in the 4 requested positions by calling the render_object function, which is located in the camera_functions.py file, after applying translation/rotation where needed.

In the second file, all the functions requested in the assignment are implemented and will be discussed in the next section of this report.

Finally, the third file is a copy of the same file that was submitted in the first assignment, with only one change. Specifically, at line 85, a conditional check using if has been added to verify if the point we intend to color is within the boundaries. If not, we choose not to color it.

## Tasks
### Function: Affine transformation
Implement the function:  

cq = affine_transform(cp, θ, u, t)

which transforms a point cp ∈ R3 (in non-homogeneous form) by rotating it by an angle θ about an axis passing through the origin of the coordinate system and parallel to u, and shifting it according to the displacement vector t, where:
* cp ∈ R3 is a 3 × 1 vector with the coordinates of a point p with respect to a coordinate system.
* θ is the rotation angle.
* u is a 3 × 1 (unit) vector representing the rotation axis.
* t is the 3 × 1 translation vector.
  
The function will implement an affine point transformation. Make sure that the affine_transform function works correctly even when cp and cq are 3 × N arrays with the coordinates of points. 

Note: The function should be able to apply either a rotation transformation or a translation transformation only, depending on the arguments (for example, if t = None, then it should apply only a rotation transformation, and similarly for translation).

### Function: Coordinate System Transformation
Implement the function:

dp = system_transform(cp, R, c0)

where:
* cp ∈ R3 is the 3 × 1 column with the coordinates of a point p with respect to a coordinate system.
* dp ∈ R3 are the coordinates of the same point in a new coordinate system with origin o ⊕ v0, which results from a rotation transformation R ∈ R3×3.
* R ∈ R3×3 is a rotation matrix.
* c0 ∈ R3 is the 3 × 1 column with the coordinates of the vector v0 with respect to the initial coordinate system.

### Function: project camera
Let cp ∈ R3 be the 3 × 1 column with the coordinates of a point with respect to the WCS {o, x0, y0, z0}. Also, suppose that a perspective camera has a center c = o ⊕ vc and unit vectors {xc, yc, zc}. Implement the function:

[verts2D, depth] = project_cam(f, cv, cx, cy, cz, p)

where cv, cx, cy, cz are the coordinates of vc, xc, yc, zc respectively, with respect to the WCS, and f is the distance from the center (measured in the units used by the camera's coordinate system).

The function will generate the perspective projections of the three-dimensional points and return them in the verts2D array with dimensions 2 × N. The function will also calculate the depth of each point before projecting it into 2 dimensions and return it in the depth array with dimensions 1 × N. Ensure that project_camera works correctly even when p is a 3 × N array with the coordinates of points.

Implement the function:

[verts2d, depth] = project_cam_lookat(f, corg, clookat, cup, verts3d)

This function will generate perspective projections and depth for the points in verts3d, similar to the previous function. However, it takes as arguments the coordinates clookat and cup (in non-homogeneous form) of the target point K and the unit up vector, respectively. The corg contains, as before, the coordinates of the camera center with respect to the WCS.

### Function: Rasterization
Implement the function:

vertsrast = rasterize(verts2d, imgh, imgw, camh, camw)

This function maps the coordinates of points from the system of a camera with a projection of dimensions camh × camw (in inches) to integer positions (pixels) of an image with dimensions imgh × imgw, generated as output by the camera during photography.

Note: The camera axis passes through the center of the rectangle with dimensions camh × camw, while the numbering of the camh × camw array of the image starts from bottom to top and from left to right, with values [0, ..., camw - 1] horizontally and [0, ..., camh - 1] vertically.

### Functon: Render Object
Implement the function:
I = render_object(verts3d, faces, vcolors, imgh, imgw, camh, camw, f, corg, clookat, cup)
where:

* I is the color image of dimensions imgh × imgw × 3. The image will contain K colored triangles.
* verts3d are the three-dimensional coordinates of the object's points of size L × 3.
* faces is the array containing indices to points in the verts3d array that make up the vertices of the triangles. The array has dimensions × 3. The i-th row of the array specifies the three vertices forming the triangle (with reference to vertices in the verts3d array, and numbering starting from 0).
* vcolors is the array with the colors of the vertices. The vcolors array has dimensions L × 3. The i-th row of the array specifies the color components of the corresponding vertex.
* imgh and imgw are the height and width of the canvas, respectively.
* camh and camw are the height and width of the camera's projection (in inches).
* f is the distance from the projection to the center (measured in the units used by the camera's coordinate system).
* clookat is the coordinates of the target point (in non-homogeneous form).
* cup is the unit up vector of the camera (in non-homogeneous form).
* corg contains the coordinates of the camera's center with respect to the WCS.
  
The function should use the above functions appropriately to implement the entire rendering pipeline of an object. It should also use the render function from the previous task to color the object using the Gouraud shading method.

## Instructions
Create the script for demonstration, named demo.py. The script should be called without external arguments, read the object from the file hw2.npy provided to you, and execute a predefined set of transformations described below:

As input, use the array verts3d, which contains the three-dimensional coordinates of the K vertices of the triangles that make up the object. Given the points of the verts3d array, your script should sequentially perform the following steps:

(a) Translate them by t1.
(b) Rotate them by an angle ϕ radians around an axis parallel to the vector u.
(c) Translate them by t2.

Each step takes the output of the previous step as input. After each step, you should capture an image of the object by calling the render_object function with camera parameters corg, clookat, cup, and color it by calling the render function from the first assignment using Gouraud shading.

Overall, you should generate 4 images of the object: one in its initial position and one for the results of steps (a) - (c). Each image should be saved with a filename corresponding to the step number (considering that the initial position is step 0) and have the .jpg extension.

Note: The file hw2.npy provided to you contains the parameters of the object (verts_3d, vcolors, faces), as well as all known parameters (camera parameters, translation vectors, rotation axes, etc.).

## Assumptions
* imgw = imgh = 512
* camw = camh = 15
* f = 70
* Do not rotate the object in order to look straight

## Implementations
### Function: Calclulate Rotation Matrix
R = calculateRotMatrix(theta,u)

This function is optional and was added for easier code understanding. It takes as input the angle in radians (theta) at which we want the rotation to occur and the axis of rotation (u). First, we check the case where no rotation is needed, in which case we return a 3x3 identity matrix. In all other cases, we calculate the rotation matrix R and return it based on the known formula from the theory:

![Rotation Matrix Calculation](/images/rot_matrix.png)

### Function: Affine transformation
First, we calculate the rotation matrix using the previous function. Then, we check if we need a translation according to the input vector t. Because we are using homogeneous coordinates, we need to construct the 4x4 matrix in every case. The function works properly even when the input c_p has dimensions of 3xN instead of 3x1.

### Function: System transformation
We follow a similar process as before, but now we need to multiply the translation vector co by the matrix -R, again following the methodology presented in the lectures of the course. Also we take care to give as input in the function the inverse matrix of R.

![Coordinate System Transformation](/images/system_transform.png)

### Function: Project Camera
We define the rotation matrix based on the 3 vectors c_x, c_y, c_z, which constitute the camera's basis vectors, and pass its inverse as input to the system_transform function. Here, we know that if a matrix is a rotation matrix, its inverse is equal to its transpose. After transforming the coordinates to the CCS (Camera Coordinate System), we calculate the matrix with the depth of each point as its dimension along the z-axis and find the (x, y) coordinates by dividing them by the depth and multiplying by the scaling factor f, as defined by the theory.

### Function: Rasterization
We calculate the coefficients k, which are used for the transformation from the coordinates in the projection to the image. Additionally, we know that the center of the projection is at coordinates cam_h/2,cam_w/2. Therefore, we add these values to the initial coordinates and then multiply by the coefficients k to obtain the final result, which is returned.

### Function: Render Object
The function that connects all the previous ones to produce the final photograph. Initially, the function project_cam_lookat is called, which will return the coordinates of the points and their depths. Then, the rasterize function will transform the coordinates from those of the projection into an image, which will be used in the render function from the previous assignment to render the complete image and return it as the output of the function.

## Results
1. Image without processing, from the original data:

![Result 0](/images/result0.jpg)

2. Image after shift=t1

![Result 1](/images/result1.jpg)

3. Image after rotation φ angles

![Result 2](/images/result2.jpg)

4. Image after shift=t2

![Result 3](/images/result3.jpg)







