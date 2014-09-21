#include <cstddef>
#include <cmath>

#include "Extra Math.h"
#include "Runge Kutta.h"
#include "Runge Kutta Data.h"
#include "Differential Geometry.h"
#include "Files.h"
#include "Camera.h"
#include "main.h"

int lightaxes;

//Construct the camera;
Camera3D::Camera3D(BlackHole *bh, double *position, double *look, int coord, int cameratype, int pixh, int pixv, double fovh, double fovv, double zoomh1, double zoomh2, double zoomv1, double zoomv2)
{
    directions[1]  = 1.0;
    directions[2]  = 0.0;
    directions[3]  = 0.0;
    directions[4]  = 0.0;
    directions[5]  = 0.0;
    directions[6]  = 1.0;
    directions[7]  = 0.0;
    directions[8]  = 0.0;
    directions[9]  = 0.0;
    directions[10] = 0.0;
    directions[11] = 1.0;
    directions[12] = 0.0;
    directions[13] = 0.0;
    directions[14] = 0.0;
    directions[15] = 0.0;
    directions[16] = 1.0;

    blackhole = bh;
    this->fovh = fovh;
    this->fovv = fovv;
    this->zoomh1 = zoomh1;
    this->zoomh2 = zoomh2;
    this->zoomv1 = zoomv1;
    this->zoomv2 = zoomv2;
    this->pixh = pixh;
    this->pixv = pixv;

    GenerateCamera(position, look, coord, cameratype);
}

Camera3D::~Camera3D()
{
    delete[] light;
    delete blackhole;
}

//Takes the black hole mass and angular momentum, the position of the camera, the look matrix,
//  an integer code specifying the coordinate types and an integer code specifying the type
//  of camera. 
//The coordinates: 
//  0 - past EF coordinates.
//  1 - future EF coordinates.
//  2 - cartesian r > 0.
//  3 - cartesian r < 0.
//The camera types: 
//  0 - Spherical camera, like the human eye.
//  1 - Flat camera, like an ordinary camera. 
void Camera3D::GenerateCamera(double *position, double *look, int coord, int cameratype)
{
    double newvec[12];
    double *cameraup, *cameraleft, *camera;
    double length, x, y, z, theta, phi;
    static int i, j,k;
    pix = pixh * pixv;
    light = new double[pix*13];

    double metric[16];
    coordinate = CoordinateTest(coord, position);
    if(coord == 3) blackhole->reverseMass();
    if(coordinate != coord) 
        CoordinateSwitch(blackhole, position, look, 4, coordinate, coord);

    switch (coordinate) {
        case PAST_EF :
            blackhole->MetricPastEddingtonFinkelstein(position, metric);
            break;
        case FUTURE_EF:
            blackhole->MetricFutureEddingtonFinkelstein(position, metric);
            break;
        case CARTESIAN_POSITIVE_R:
            blackhole->MetricCartesian(position, metric);
            break;
        case CARTESIAN_NEGATIVE_R:
            blackhole->reverseMass();
            blackhole->MetricCartesian(position, metric);
            break;
    }

    for(i = 0; i < 16; i++) directions[i] = look[i];
    pos[0] = position[0]; pos[1] = position[1]; pos[2] = position[2]; pos[3] = position[3]; 
    Orthonormalize(directions, metric);

    Transpose(directions);

    for(i = 0; i < pixh; i++) {
        for(j = 0; j < pixv; j++) {
            camera = &light[(pixv*i+j)*13];
            cameraup = &camera[4];
            cameraleft = &camera[8];

            switch (cameratype) {
                case CAMERA_SPHERICAL:
                    theta = 0.785398*(-fovv + 2*fovv * (zoomv1 + (zoomv2 - zoomv1) * (j + .5) /((double) pixv)));
                    phi = 0.785398*(-fovh + 2*fovh * (zoomh1 + (zoomh2 - zoomh1) * (i + .5) /((double) pixh)));
                    z = cos(theta)*cos(phi);
                    x = -cos(theta)*sin(phi);
                    y = sin(theta);
                    break;
                case CAMERA_FLAT:
                    camera[12] = 1.0;
                    x = -fovh + 2*fovh * (zoomh1 + (zoomh2 - zoomh1) * (i + .5) /((double) pixh));
                    y = -fovv + 2*fovv * (zoomv1 + (zoomv2 - zoomv1) * (j + .5) /((double) pixv));
                    z = 1.0;
            }

            length = pow(x*x + y*y + z*z, .5);
            newvec[3] = x/length; newvec[1] = z/length; newvec[2] = y/length;
            newvec[0] = -pow(newvec[1]*newvec[1] + newvec[2]*newvec[2] + newvec[3]*newvec[3], .5);
            MatrixDotVector(camera, directions, newvec);

            length = pow(y*y + 1, .5);
            newvec[3] = 0;
            newvec[1] = -y/length;
            newvec[2] = 1/length;
            newvec[0] = -pow(newvec[1]*newvec[1] + newvec[2]*newvec[2] + newvec[3]*newvec[3], .5);
            MatrixDotVector(cameraup, directions, newvec);

            length = pow((1 + x*x + y*y)/(y*y + 1), .5);
            newvec[3] = -1/length;
            newvec[1] = x/((1+y*y)*length);
            newvec[2] = y*newvec[1];
            newvec[0] = -pow(newvec[1]*newvec[1] + newvec[2]*newvec[2] + newvec[3]*newvec[3], .5);
            MatrixDotVector(cameraleft, directions, newvec);

        }
    }
}

//Takes a snapshot with the given camera. 
//Arguments: 
//  The Runge-Kutta method to be used.
//  The filename.
//  Maximum number of steps a light ray may take.
//  The maximum distance from the black hole a ray is allowed to go before it is stopped. 
//  The initial step size.
//  Coordinate code specifying the output coordinates to say where the light ray ends up. 
//  A coordinate test function. 
void Camera3D::Snapshot(RungeKuttaForm const *solvedata, char const*file, int maxsteps, double maxdistance, double initstep, double errscale, int outputcoordinate, int (*coordtest)(int, double[], double[]))
{
    double error[17];
    double *lambda1, *lambda2, *temp, *dlambda1, *dlambda2;
    lambda1 = new double[17]; lambda2 = new double[17]; dlambda1 = new double[17]; dlambda2 = new double[17];
    double step, s, dist, distmin, distmax, previousdist, *v;
    int presentcoordinate, newcoordinate;
    ofstream out;
    out.open(file);
    int i, j, counter = 0, outeradjustments;
    double maxdistancesquared = maxdistance*maxdistance;
    double stepmax, stepmin;
    bool initeval;
    double inversedirections[16];
    void (*initialderivatives)(double, double[], double[]);
    void (*derivatives)(double, double[], double[]);
    double (*initialdistance)(double[]);
    double (*distance)(double[]);
    RungeKuttaData *solver = new RungeKuttaData(solvedata, 8);

    //The different coordinate schemes:
    //0: past EF coordinates
    //1: future EF coordinates
    //2: cartesian r > 0.
    //3: cartesian r < 0.
    switch (coordinate){
        case PAST_EF:
            initialderivatives = DerivativesPastEddingtonFinkelstein;
            initialdistance = DistancePastEddingtonFinkelstein;
            break;

        case FUTURE_EF:
            initialderivatives = DerivativesFutureEddingtonFinkelstein;
            initialdistance = DistanceFutureEddingtonFinkelstein;
            break;

        case CARTESIAN_POSITIVE_R:
            initialderivatives = DerivativesCartesian;
            initialdistance = DistanceCartesian;
            break;

        case CARTESIAN_NEGATIVE_R:
            initialderivatives = DerivativesCartesian;
            initialdistance = DistanceCartesian;
            break;
    }

    InverseMatrix(inversedirections, directions);

    out << pixh << " " << pixv << endl;

    for(i=0;i<pix;i++)
    {
        s = 0.0; step = initstep;
        dist = 0.0;
        for(j = 0; j < 4; j++) lambda1[j] = pos[j];
        while(j < 8) {lambda1[j] = light[13*i+j-4];j++;}
        for(j = 0; j < 8; j++) lambda2[j] = lambda1[j];

        derivatives = initialderivatives;
        distance = initialdistance;
        initeval = true;
        presentcoordinate = coordinate;

        for(j=0;j<maxsteps;j++)
        {
            counter ++;

            if(initeval)derivatives(0.0, lambda1, dlambda1);
            initeval = false;

            step = solver->timestep(derivatives, lambda1, dlambda1, 0.0, step, 1., errscale, lambda2, error, dlambda2);
            if (((step >= -.0000001)&&(step <= .0000001)) || isnan(step) || isinf(step))   break;

            distmin = dist;
            dist = distance(lambda2);

            //If the light ray goes outside the max distance, it is necessary to track it
            //  to find its location at exactly that distance. This loop tests different
            //  lengths for the last time step until the light ray is as close as possible
            //  to the max distance. IF this is not done, the picture comes out a little 
            //  funny-looking. 
            if (dist > maxdistancesquared)
            {
                outeradjustments = 0;
                previousdist = pow(dist, .5);
                distmin = pow(distmin, .5);
                distmax = previousdist;
                stepmax = step; stepmin = 0;
                while(true)
                {
                    dist = pow(dist, .5);

                    if(dist < maxdistance) {
                        stepmin = step;
                        distmin = dist;
                    } else {
                        if(dist - DISTANCERANGE < maxdistance) break;
                        else {
                            stepmax = step; distmax = dist;
                        }
                    }

                    step = stepmin + (stepmax - stepmin)*(maxdistance - distmin)/(distmax - distmin);
                    solver->step(derivatives, lambda1, dlambda2, 0.0, step, lambda2, error, dlambda2);

                    dist = distance(lambda2);

                    if(outeradjustments > MAXOUTERADJUSTMENTS) break; 
                    outeradjustments++;
                }
                break;
            }

            v = &lambda1[4];
            //Test if a coordinate change is needed. 
            newcoordinate = coordtest(presentcoordinate, lambda2, lambda2);
            if(newcoordinate != presentcoordinate) {
                switch (newcoordinate){
                    case PAST_EF:
                        derivatives = DerivativesPastEddingtonFinkelstein;
                        distance = DistancePastEddingtonFinkelstein;
                        break;
                    case FUTURE_EF:
                        derivatives = DerivativesFutureEddingtonFinkelstein;
                        distance = DistanceFutureEddingtonFinkelstein;
                        break;
                    case CARTESIAN_POSITIVE_R:
                        derivatives = DerivativesCartesian;
                        distance = DistanceCartesian;
                        break;
                    case CARTESIAN_NEGATIVE_R:
                        derivatives = DerivativesCartesian;
                        distance = DistanceCartesian;
                        break;
                }
                CoordinateSwitch(blackhole, lambda2, v, 1, newcoordinate, presentcoordinate);
                initeval = true;
            }
            presentcoordinate = newcoordinate;

            temp = lambda1; lambda1 = lambda2; lambda2 = temp;
            temp = dlambda1; dlambda1 = dlambda2; dlambda2 = temp;
        }
        CoordinateSwitch(blackhole, lambda2, &lambda2[4], 1, 4, presentcoordinate);

    }
    delete solver;
    delete[] lambda1; delete[] lambda2; delete[] dlambda1; delete[] dlambda2;
}
