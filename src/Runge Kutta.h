#pragma once

//The RungeKuttaForm contains the abstract data about a runge-kutta method.
//It cannot be used to solve a differential equation. 
class RungeKuttaForm{
public:
    int steps; //The number of steps.
    int order; //The order of the method.
    bool fsal; //First-same-as-last property.
    double *a; //The data for the factors on the intermediate steps.
    double *b; //The data for the factors on the final step.
    double *c; //The data for the factors on the time steps.
    double *d; //The data for the factors on the error estimation.
    RungeKuttaForm(int st, int o, bool f, double rka[], double rkb[], double rkc[], double rkd[]);
    ~RungeKuttaForm();
};

//The RungeKuttaData knows how to run a differential equation and includes
//an extra parameter 'n', which determines the level of precision. 
class RungeKuttaData{
    int steps;   //The number of steps.
    int order;   //The order of the method.
    int n;       //The length of the arrays.
    bool fsal;   //First-same-as-last property.
    double *a;   //The data for the factors on the intermediate steps.
    double *b;   //The data for the factors on the final step.
    double *c;   //The data for the factors on the time steps.
    double *db;  //The data for the factors on the error estimation.
    double *K;   //An array for storing temporary data during the computation.
    double *tmp; //Another temporary storage value.
public:
    RungeKuttaData(RungeKuttaForm *data, int nn);
    ~RungeKuttaData();

    void step(void(*f)(double, double [], double []), double x[], double dx[], double t, double dt, double x2[], double err[], double dx2[]);

    double timestep(void (*f)(double, double [], double[]), double x[], double a[], double t, double dt, double eps, double errscale, double x2[], double err[], double dx2[]);
};

//double rkstepsizer(void (*f)(double, double [], double[]), RungeKuttaData *data, double x[], double a[], int n, double t, double dt, double eps, double errscale, double x2[], double err[], double dx2[]);
