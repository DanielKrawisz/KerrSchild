#pragma once

#define PAST_EF              0
#define FUTURE_EF            1
#define CARTESIAN_POSITIVE_R 2
#define CARTESIAN_NEGATIVE_R 3
#define RADIAL               4

class BlackHole
{
private:
    double irad; //the inner horizon.
    double orad; //the outer horizon.

    double M;     //Black hole mass
    double a;     //The rotational parameter
    bool horizon; //Whether there is a horizon.
public:

    BlackHole(double bh_M, double bh_a);
    ~BlackHole();

    double mass() {return M;}
    double rotation() {return a;}
    bool horizon_exists() {return horizon;}
    void reverseMass() {M = -M;}

    //These functions set g as the black hole metric tensor for
    //  various kinds of coordinates at a given location. 
    void MetricPastEddingtonFinkelstein(double x[], double g[]);
    void MetricFutureEddingtonFinkelstein(double x[], double g[]);
    void MetricCartesian(double x[], double g[]);

    //These set dx as the acceleration of an object at the location
    //  x for various kinds of black hole coordinates. 
    void AccelerationPastEddingtonFinkelstein(double x[], double dx[]);
    void AccelerationFutureEddingtonFinkelstein(double x[], double dx[]);
    void AccelerationCartesian(double x[], double dx[]);
};

//These are the functions passed to the numerical differential
//  equation solver. The first argument is unused. The second is 
//  the coordinates of a light ray. The third is set to the acceleration
//  vector of the light ray at that point. 
void DerivativesPastEddingtonFinkelstein(double t, double x[], double dx[]);
void DerivativesFutureEddingtonFinkelstein(double t, double x[], double dx[]);
void DerivativesCartesian(double t, double x[], double dx[]);

//The distance from the black hole in various kinds of coordinates.
double DistancePastEddingtonFinkelstein(double x[]);
double DistanceFutureEddingtonFinkelstein(double x[]);
double DistanceCartesian(double x[]);

//Different kinds of coordinates are convenient in different regions of 
//  the space around a black hole. The differential equation solver
//  requires a test to determine when the light ray has entered a
//  region in which a new kind of coordinate is convenient. 
//This is what CoordinateTest does. It takes an integer representing
//  the coordinates already being used and returns an integer
//  represeting the coordinates that should be used next.
//The integer codes are: 
//  0 - past EF coordinates.
//  1 - future EF coordinates.
//  2 - cartesian r > 0.
//  3 - cartesian r < 0.
//The version that takes two vectors is used during the simulation when
//  the direction of light determines which coordinate should be used. 
//The version that only takes one vector is for the very start of the
//  simulation to determine what coordinates the initial light ray should
//  be represented in. 
int CoordinateTest(int coord, double *prev, double *next);
int CoordinateTest(int coord, double *prev);

//Given the integer code returned by one of the CoordinateTest functions, the 
//  CoordinateSwitch function implements the necessary change of coordinates.
void CoordinateSwitch(BlackHole *blackhole, double lambda2[], double v[], int n, int newcoordinate, int presentcoordinate);
