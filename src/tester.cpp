#ifndef TESTER_CPP
#define TESTER_CPP

#include "main.h"
using namespace std;
using namespace RungeKutta;

void Tester(){
    //RungeKuttaTester();
    //AccelerationTester();
    CoordinateTester();
}

void AccelerationTester(){
    double x[8] = {0.0, 7.0, 8.0, 9.0, 1.0, 0.0, 0.0, 0.0};
    double dx[8];
    int i;
    BlackHole *blackhole = new BlackHole(2., 1.9);

    cout << "Testing the acceleration of the Cartesian coordinate frame: " << endl;
    for(i = 0; i < 8; i++){
        cout << x[i] << " ";
    } cout << endl;
    DerivativesCartesian(0, x, dx);
    cout << "Has acceleration " << endl;
    for(i = 0; i < 8; i++){
        cout << dx[i] << " ";
    } cout << endl;
    cout << "x: " << x[0] << " " << x[1] << endl;
    cout << "Now switching to past Eddington-Finkelstein." << endl;
    //CartesianToRadialCoordinateTransformation(&x[4], x, blackhole->M, blackhole->a, 1);
    CartesianToRadialCoordinateConversion(x, blackhole->M, blackhole->a, 1);
    EddingtonFinkelsteinCoordinateTransformation(&x[4], x, blackhole->M, blackhole->a, -1);
    EddingtonFinkelsteinCoordinateConversion(x, blackhole->M, blackhole->a, -1);
    for(i = 0; i < 8; i++){
        cout << x[i] << " ";
    } cout << endl;
    DerivativesPastEddingtonFinkelstein(0, x, dx);
    cout << "Has acceleration " << endl;
    for(i = 0; i < 8; i++){
        cout << dx[i] << " ";
    } cout << endl;

    //delete blackhole;
}

void RungeKuttaTester(){
    double x[2] = {1.0, 0.0}, x_[2], v[2], v_[2], err[2];
    RungeKuttaForm *RKdata[] = {Euler, Midpoint, Kutta, BogackiSchampine, Fehlberg, CashKarp, DormandPrince};
    RungeKuttaData *RK;
    cout << setiosflags(ios::left);
    cout << "Runge Kutta Test." << endl;
    for(int i = 0; i < 7; i++) {
        cout << IND << "Test " << i << endl;
        QuarticOscillator(0.0, x, v);
        cout << IND << IND << setw(ITEM) << "input:" << setw(ITEM) << x[0] << setw(ITEM) << x[1] <<  endl;
        cout << IND << IND << setw(ITEM) << "velocity:" << setw(ITEM) << v[0] << setw(ITEM) << v[1] << endl;
        RK = new RungeKuttaData(RKdata[i], 2);
        RK->step(QuarticOscillator, x, v, 0.0, 0.1, x_, err, v_);
        cout << IND << IND << setw(ITEM) << "next:" << setw(ITEM) << x_[0] << setw(ITEM) << x_[1] << endl;
        cout << IND << IND << setw(ITEM) << "error:" << setw(ITEM) << err[0] << setw(ITEM) << err[1] << endl;
        delete RK;
    }
    cout << "Time step test." << endl;
    for(int i = 0; i < 7; i++){
        double time;
        cout << IND << "Test " << i << endl;
        QuarticOscillator(0.0, x, v);
        cout << IND << IND << setw(ITEM) << "input:" << setw(ITEM) << x[0] << setw(ITEM) << x[1] <<  endl;
        cout << IND << IND << setw(ITEM) << "velocity:" << setw(ITEM) << v[0] << setw(ITEM) << v[1] << endl;
        RK = new RungeKuttaData(RKdata[i], 2);
        time = RK->timestep(QuarticOscillator, x, v, 0.0, 0.1, .001, 1.0, x_, err, v_);
        cout << IND << IND << setw(ITEM) << "next:" << setw(ITEM) << x_[0] << setw(ITEM) << x_[1] << endl;
        cout << IND << IND << setw(ITEM) << "step:" << setw(ITEM) << time << endl;
        delete RK;
    }
}

void CoordinateTester(){
    double M = 2.0;
    double a = 1.9;
    double ovecs[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    double vecs[4][4];
    double pos[][4] =  {{0, 7, 8, 9}, {0, 2, 0, .1}, {0, .3, .4, .5}, {0, 7, 8, 9}};
    int i, j, k;
    cout << "Coordinate tester." << endl;

    struct carp{
        static void printvectors(int i, double vecs[][4], double pos[][4]){
            cout << IND;
            int j; int k;
            for(j = 0; j < 4; j++) cout << setw(ITEM) << pos[i][j];
            cout << endl;
            for(k = 0; k < 4; k++){
                cout << IND;
                for(j = 0; j < 4; j++) cout << setw(ITEM) << vecs[k][j];
                cout << endl;
            }
        }
    };

    for(i = 0; i < 4; i++){
        for(k = 0; k < 4; k++){
            for(j = 0; j < 4; j++){
                vecs[k][j] = ovecs[k][j];}}
        cout << "Location " << i << endl;
        cout << "original: " << endl;
        carp::printvectors(i, vecs, pos);
        cout << "Shift to radial." << endl;
        for(j = 0; j < 4; j++)
            CartesianToRadialCoordinateTransformation(&vecs[j][0],&pos[i][0], M, a, i==3?-1:1);
        CartesianToRadialCoordinateConversion(&pos[i][0], M, a, i==3?-1:1);
        carp::printvectors(i, vecs, pos);
        cout << "Shift to Eddington Finkelstein."<< endl;
        for(k = 0; k < 4; k++){
            for(j = 0; j < 4; j++){
                vecs[k][j] = ovecs[k][j];}}
        for(j = 0; j < 4; j++)
            EddingtonFinkelsteinCoordinateTransformation(&vecs[j][0],&pos[i][0], M, a, -1);
        EddingtonFinkelsteinCoordinateConversion(&pos[i][0], M, a, -1);
        carp::printvectors(i, vecs, pos);
        cout << "Shift back to radial." << endl;
        for(k = 0; k < 4; k++){
            for(j = 0; j < 4; j++){
                vecs[k][j] = ovecs[k][j];}}
        for(j = 0; j < 4; j++)
            EddingtonFinkelsteinCoordinateTransformation(&vecs[j][0],&pos[i][0], M, a, 1);
        EddingtonFinkelsteinCoordinateConversion(&pos[i][0], M, a, 1);
        carp::printvectors(i, vecs, pos);
        cout << "Shift back to Cartesian." << endl;
        for(k = 0; k < 4; k++){
            for(j = 0; j < 4; j++){
                vecs[k][j] = ovecs[k][j];}}
        for(j = 0; j < 4; j++)
            RadialToCartesianCoordinateTransformation(&vecs[j][0],&pos[i][0], M, a);
        RadialToCartesianCoordinateConversion(&pos[i][0], M, a);
        carp::printvectors(i, vecs, pos);
    }
}

void QuarticOscillator(double t, double x[], double dx[]){
    static double m = 1.0, k = 1.0, q = 0.1;
    static double x2;
    x2 = x[0]*x[0];
    dx[0] = x[1];
    dx[1] = -k*x2 - q*x2*x2;
}

#endif