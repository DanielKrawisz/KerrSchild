#include <cmath>

#include "Coordinate Transformations.h"
#include "Differential Geometry.h"
#include "Extra Math.h"

using namespace std;

//Two special functions used several times in coordinate transformations.
void BlackHoleArcTans(double &L1, double &L2, double r, double m, double a){
    static double L, S, J, Q;
    L = m*m - a*a;
    if(L > 0){
        S = sqrt(L);
        J = (r - m) / S;
        Q = (log(abs(1 - J)) - log(abs(1 + J)))/(2 * S);
    } else if(L < 0) {
        S = sqrt(-L);
        Q = atan((r - m) / S) / S;
    } else {
        Q = 1/(m - r);
    }
    L1 = a * Q; L2 = 2*m*m*Q + r + m*log(a*a - 2*m*r + r*r);
}

double BlackHoleArcTan(double r, double m, double a){
    static double L, S, J;
    L = m*m - a*a;
    if(L > 0){
        S = sqrt(L);
        J = (r - m) / S;
        return a *(log(abs(1 - J)) - log(abs(1 + J)))/(2 * S);
    } else if(L < 0) {
        S = sqrt(-L);
        return a*atan((r - m) / S) / S;
    } else {
        return m/(m - r);
    }
}

//Convert a coordinate point between kinds of Eddington-Finkelstein coordinates. The parameter 'n'
//  can be used to convert between different kinds of EF coordinates.
//  For example, a value of -1 converts back from EF coordinates to
//  ordinary radial coordinates. 
void EddingtonFinkelsteinCoordinateConversion(double *x, double m, double a, int n)
{
    static double L1, L2;

    BlackHoleArcTans(L1, L2, x[1], m, a);

    x[0] += n*L2;
    x[3] += n*L1;
}

//Convert a vector at a point between different kinds of Eddington-Finkelstein coordinates.
//  The parameter 'n' can be used to convert between different kinds of EF coordinates.
//  For example, a value of -1 converts back from EF coordinates to
//  ordinary radial coordinates. 
void EddingtonFinkelsteinCoordinateTransformation(double *v, double *x, double m, double a, int n)
{
    static double delta;

    delta = a*a - 2*m*x[1] + x[1]*x[1];

    v[0] += n*(a*a + x[1]*x[1])*v[1] / delta;
    v[3] += n*a*v[1] / delta;
}

//Convert a coordinate point from Radial to Cartesian.
void RadialToCartesianCoordinateConversion(double *x, double m, double a)
{
    static double t, r, theta, phi;
    static double Sr, St;
    static double sintheta;
    static double phifunction;
    static double sinphif, cosphif;

    t=x[0]; r=x[1]; theta=x[2]; phi=x[3];
    BlackHoleArcTans(Sr, St, r, m, a);
    sintheta = sin(theta);
    phifunction = phi + Sr;
    sinphif = sin(phifunction); cosphif = cos(phifunction);
    x[0] = t + St - r;
    x[1] = sintheta*(r*cosphif - a*sinphif);
    x[2] = sintheta*(a*cosphif + r*sinphif);
    x[3] = r*cos(theta);
}

//Convert a vector at a point from Radial to Cartesian.
void RadialToCartesianCoordinateTransformation(double *v, double *x, double m, double a)
{
    static double t, r, theta, phi, S;
    static double transformation[] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    static double delta, Sd;
    static double sintheta, costheta;
    static double phifunction;
    static double sinphif, cosphif;
    t=x[0]; r=x[1]; theta=x[2]; phi=x[3];
    S=BlackHoleArcTan(r, m, a);
    delta = a*a - 2*m*r + r*r;
    Sd = a / delta;
    sintheta = sin(theta); costheta=cos(theta);
    phifunction = phi + S;
    sinphif = sin(phifunction); cosphif = cos(phifunction);

    transformation[1 ] = (r*r + a*a)/delta - 1;
    transformation[5 ] = sintheta*((1 - a*Sd)*cosphif - r*Sd*sinphif);
    transformation[6 ] = costheta*(r*cosphif - a*sinphif);
    transformation[7 ] = sintheta*(-r*sinphif - a*cosphif);
    transformation[9 ] = sintheta*((1 - a*Sd)*sinphif + r*Sd*cosphif);
    transformation[10] = costheta*(r*sinphif + a*cosphif);
    transformation[11] = sintheta*(r*cosphif - a*sinphif);
    transformation[13] = costheta;
    transformation[14] = -r*sintheta;

    MatrixDotVector(v, transformation, v);
}

//Cartesian to Radial.
void CartesianToRadialCoordinateConversion(double *x, double m, double a, int direction)
{
    static double ra, rasq, r, theta, sintheta, rraa, sine, phi, sqr2, Sr, St;
    ra = x[1]*x[1]+x[2]*x[2]+x[3]*x[3]-a*a;
    rasq = sqrt(4*a*a*x[3]*x[3]+ra*ra);
    sqr2 = sqrt(2.);
    r = direction*sqrt((ra + rasq))/sqr2;
    BlackHoleArcTans(Sr,St,r,m,a);
    theta = acos(x[3]/r);
    sintheta = sin(theta); rraa = a*a+r*r;
    sine = (-a*x[1]+r*x[2])/(sintheta*rraa);
    phi = acos((r*x[1]+a*x[2])/(sintheta*rraa));
    if(sine < 0) {
        phi *=-1;
        //dphi *=-1;
    }

    x[0] += r - St;
    x[1] = r;
    x[2] = theta;
    x[3] = phi - Sr;
}

//Convert a vector at a point from Cartesian to Radial.
void CartesianToRadialCoordinateTransformation(double *v, double *x, double m, double a, int direction)
{
    //TODO: figure out all the cases for what to do with the logorithms with negative and positive r. 
    static double conversion[] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    static double roa, roz, r, rdx, rdy, rdz, rq, th, thdx, thdy, thdz, r2sq, dtr, sign,
                    a2r2, a2r22, sinth, rxay, delta, juju, jujuju, jujujuju, dphr, dphth;
    roa = x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - a*a;
    roz = sqrt(4*a*a*x[3]*x[3] + roa*roa);
    r = direction*sqrt((roa + roz)/2); rq = r/roz;
    rdx = x[1]*rq; rdy = x[2]*rq; rdz = x[3]*(1 + (roa + 2*a*a)/roz)/(2*r);
    r2sq = r*r*sqrt(1 - x[3]*x[3]/(r*r));
    th = acos(x[3]/r); thdx = x[3]*rdx/r2sq; thdy = x[3]*rdy/r2sq; thdz = (rdz*x[3] - r)/r2sq;
    a2r2 = a*a + r*r; a2r22 = a2r2*a2r2; 
    delta = a2r2 -2*m*r;
    dtr = 1-a2r2/delta;
    sign = (-a*x[1] + r*x[2])?1:-1;
    sinth = sin(th); rxay = r*x[1] + a*x[2]; 
    juju = rxay*rxay/(a2r22*sinth*sinth);
    jujuju = sqrt(1-juju);
    jujujuju = jujuju*a2r2*sinth;
    dphr = -a/delta - sign*(x[1]/(sinth*a2r2) - 2*r*rxay/(sinth*a2r22))/jujuju;
    dphth = sign*rxay*cot(th)/jujujuju;
    conversion[1 ] = dtr*rdx;
    conversion[2 ] = dtr*rdy;
    conversion[3 ] = dtr*rdz;
    conversion[5 ] = rdx;
    conversion[6 ] = rdy;
    conversion[7 ] = rdz;
    conversion[9 ] = thdx;
    conversion[10] = thdy;
    conversion[11] = thdz;
    conversion[13] = dphr*rdx + dphth*thdx - sign*r/jujujuju;
    conversion[14] = dphr*rdy + dphth*thdy - sign*a/jujujuju;
    conversion[15] = dphr*rdz + dphth*thdz;

    MatrixDotVector(v, conversion, v);

}
