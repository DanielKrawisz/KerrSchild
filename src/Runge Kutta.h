#pragma once

//The RungeKuttaForm contains the abstract data about a runge-kutta method.
//It cannot be used to solve a differential equation. 
class RungeKuttaForm{
private: 
    int steps; //The number of steps.
    int order; //The order of the method.
    bool fsal; //First-same-as-last property.
    double *a; //The data for the factors on the intermediate steps.
    double *b; //The data for the factors on the final step.
    double *c; //The data for the factors on the time steps.
    double *d; //The data for the factors on the error estimation.
public:
    int get_steps() const {return steps;}
    int get_order() const {return order;}
    bool first_same_as_last() const {return fsal;}
    double const *get_parameter_a() const {return a;}
    double const *get_parameter_b() const {return b;}
    double const *get_parameter_c() const {return c;}
    double const *get_parameter_d() const {return d;}
    double get_parameter_a(int i) const {return a[i];}
    double get_parameter_b(int i) const {return b[i];}
    double get_parameter_c(int i) const {return c[i];}
    double get_parameter_d(int i) const {return d[i];}
    RungeKuttaForm(int st, int o, bool f, double rka[], double rkb[], double rkc[], double rkd[]);
    ~RungeKuttaForm();
};

//The RungeKuttaData knows how to run a differential equation and includes
//an extra parameter 'n', which determines the level of precision. 
class RungeKuttaData{
private:
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
    RungeKuttaData(RungeKuttaForm const *data, int nn);
    ~RungeKuttaData();

    void const step(void(*f)(double, double [], double []), double x[], double dx[], double t, double dt, double x2[], double err[], double dx2[]);

    double const timestep(void (*f)(double, double [], double[]), double x[], double a[], double t, double dt, double eps, double errscale, double x2[], double err[], double dx2[]);
};

//double rkstepsizer(void (*f)(double, double [], double[]), RungeKuttaData *data, double x[], double a[], int n, double t, double dt, double eps, double errscale, double x2[], double err[], double dx2[]);
