#pragma once

//Two special functions used several times these coordinate transformations. 
//  Not useful anywhere else. 
void BlackHoleArcTans(double &L1, double &L2, double r, double m, double a);
double BlackHoleArcTan(double r, double m, double a);

//The functions named 'transformation' take a location 'x' in some old
//  coordinates and a vector 'v' at that location given in old coordinates
//  and transform the vector so that it is in new coordinates. 
//The functions that are named 'conversion' take a location and 
//  reset it to the coordinates for that same location in
//  different coordinates. 
//Remember to use the transformation function first and THEN the 
//  conversion function! 
//The 'm' and 'a' parameters refer to the mass and rotation of the black hole
//  in whose coordinates we are transforming. 

//From radial to Eddington-Finkelstein coordinates. The parameter 'n'
//  can be used to convert between different kinds of EF coordinates.
//  For example, a value of -1 converts back from EF coordinates to
//  ordinary radial coordinates. 
void EddingtonFinkelsteinCoordinateTransformation(double *v, double *x, double m, double a, int n);
void EddingtonFinkelsteinCoordinateConversion(double *x, double m, double a, int n);

//Radial to Cartesian.
void RadialToCartesianCoordinateTransformation(double *v, double *x, double m, double a);
void RadialToCartesianCoordinateConversion(double *x, double m, double a);

//Cartesian to Radial. 
void CartesianToRadialCoordinateTransformation(double *v, double *x, double m, double a, int direction);
void CartesianToRadialCoordinateConversion(double *x, double m, double a, int direction);
