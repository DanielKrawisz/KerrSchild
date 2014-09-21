#pragma once
#include "Black Hole.h"
#include "Runge Kutta.h"

#define CAMERA_SPHERICAL            0
#define CAMERA_FLAT                 1

//The camera contains data describing the location and parameters of the camera,
//  and a link to a black hole object, which contains the data about the black
//  hole in this universe. 
class Camera3D
{
private:
    //Coordinates of the camera's location. 
    double pos[4];
    //The type of coordinates: 
    //  0 - past EF coordinates.
    //  1 - future EF coordinates.
    //  2 - cartesian r > 0.
    //  3 - cartesian r < 0.
    //  4 - radial coordinates. (standard Schwarzchild coordinates).
    int coordinate;
    //Camera parameters.
    // fovh, fovv - field of view values.
    // zoom - zoom in on any rectangle defined by these four values.
    double fovh, fovv, zoomh1, zoomh2, zoomv1, zoomv2;
    //The pixels along the horizontal and vertical, and the total number.
    int    pixh, pixv, pix;
    //A matrix defining the parameters of the camera for computational purposes.
    double directions[16];
    //The one black hole in the universe. 
    BlackHole *blackhole;
    //Stores the data describing each light ray in the camera. 
    double *light;
  
    void GenerateCamera(double *position, double *look, int coord, int cameratype);

public:
    //Takes the black hole mass and angular momentum, the position of the camera, the look matrix,
    //  an integer code specifying the coordinate types and an integer code specifying the type
    //  of camera. 
    //The coordinates: 
    //The camera types: 
    //  0 - Spherical camera, like the human eye.
    //  1 - Flat camera, like an ordinary camera. 
    Camera3D(BlackHole *bh, double *position, double *look, int coord, int cameratype, int pixh, int pixv, double fovh, double fovv, double zoomh1, double zoomh2, double zoomv1, double zoomv2);
    ~Camera3D();

    //Takes a snapshot with the given camera. 
    //Arguments: 
    //  The Runge-Kutta method to be used.
    //  The filename.
    //  Maximum number of steps a light ray may take.
    //  The maximum distance from the black hole a ray is allowed to go before it is stopped. 
    //  The initial step size.
    //  Coordinate code specifying the output coordinates to say where the light ray ends up. 
    //  A coordinate test function. 
    void Snapshot(RungeKuttaForm const* solvedata, char const* file, int maxsteps, double maxdistance, double initstep, double errscale, int outputcoordinate, int (*coordtest)(int, double[], double[]));
};
