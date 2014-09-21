#include "Runge Kutta.h"
#include "Extra Math.h"

using namespace std;

RungeKuttaForm::RungeKuttaForm(int st, int o, bool f, double rka[], double rkb[], double rkc[], double rkd[]){
//    int i;
    steps = st;
    order = o;
    fsal = f;
    a = rka; b = rkb; c=rkc; d=rkd;
}

RungeKuttaForm::~RungeKuttaForm(){
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
}

RungeKuttaData::RungeKuttaData(RungeKuttaForm const* data, int nn){
    int i;
    steps = data->get_steps()-1;
    order = data->get_order();
    fsal = data->first_same_as_last();
    n = nn;
    a = new double[steps*steps];
    b = new double[steps + 1];
    c = new double[steps];
    db = new double[steps + 1];
    K = new double[steps*n];
    tmp = new double[n];

    for(i = 0; i < steps*steps; i++) a[i] = data->get_parameter_a(i);
    for(i = 0; i < steps; i++) c[i] = data->get_parameter_c(i);
    for(i = 0; i <= steps; i++) {
        b[i] = data->get_parameter_b(i);
        db[i] = data->get_parameter_b(i) - data->get_parameter_d(i);
    }
}

RungeKuttaData::~RungeKuttaData(){
    delete[] a; delete[] b; delete[] c; delete[] db;
    delete[] K;
    delete[] tmp;
}

//One Runga-Kutta step. 
// f        - A function that takes a time, a position, and a velocity, and gives an acceleration.
// x        - the position.
// dx       - the velocity. 
// t        - the time.
// dt       - the time step.
// x2       - output: the new position.
// err      - output: The estimated error. 
// dx2      - output: the new velocity.
void const RungeKuttaData::step(void(*f)(double, double [], double []), double x[], double dx[], double t, double dt, double x2[], double err[], double dx2[]) {
    static int i, j, k, apos, kk;
    //Note that k is a counter whereas K stores data.
    apos = 0;
    //First the intermediate steps are calculated.
    for(i=0;i<steps;i++) {
        for(k=0;k<n;k++) {
            tmp[k] = 0; 
            kk = 0;
            for(j=1;j<=i;j++) {
                tmp[k] += a[apos+j]*K[kk + k];
                kk+=n;
            }
            tmp[k] = x[k] + dt*(tmp[k] + a[apos]*dx[k]);
        }
        f(t + dt*c[i], tmp, &K[n*i]);
        apos += steps;
    }
    //Then the output and error are calculated. 
    for(k=0;k<n;k++){
        x2[k] = 0;
        err[k] = 0;
        kk = 0;
        for(j=1;j<=steps;j++){
            x2[k]   += b[j]*K[kk + k];
            err[k] += db[j]*K[kk + k];
            kk+=n;
        }
        x2[k] = x[k] + dt*(x2[k] + b[0]*dx[k]);
        err[k] = dt*(err[k] + db[0]*dx[k]);
    }
    //calculate the derivative of the next step. 
    if(fsal){
        apos = n*(steps - 1);
        for(i=0;i<n;i++)
            dx2[i] = K[apos + i]; 
    } else f(t+dt, x2, dx2);

}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRINK -0.25
#define ERRCON 1.89

//From Numerical Recipies in C with some slight changes.
// f        - A function that takes a time, a position, and a velocity, and gives an acceleration.
// x        - the position.
// a        - the previous acceleration. 
// t        - the time.
// dt       - the time step.
// eps      - not used. 
// errscale - a parameter which adjusts acceptible error scales.
// x2       - output: new position.
// err      - output: estimated error.
// dx2      - output: new velocity.
double const RungeKuttaData::timestep(void (*f)(double, double [], double[]), double x[], double a[], double t, double dt, double eps, double errscale, double x2[], double err[], double dx2[])
{
    static int i;
    static double errmax, accelmax, h, htemp, xnew, *ytemp;
    h=dt; double dttemp;

    while(true)
    {
        step(f, x, a, t, dt, x2, err, dx2);
        errmax=0.0; accelmax=0.0;
        //TODO: add SAFETY and PSHRINK and PGROW. 
        for(i=0;i<n;i++) accelmax = dmax(accelmax, fabs(a[i]));
        for(i=0;i<n;i++) errmax=dmax(errmax, fabs(err[i]/(dt*errscale*accelmax)));
        if (errmax <= 1.0) break;
        dttemp=SAFETY*dt*pow(errmax,PSHRINK);
        dt=(dt >= 0.0 ? dmax(dttemp, 0.1*dt) : dmin(dttemp, 0.1*dt));

//TODO: Generate an exception. 
//		if(t + dt==t) MessageBox(0, "stepsize underflow in rkqs", "Underflow!", MB_ICONEXCLAMATION|MB_OK);
        if(t + dt == t) {/*cout << "stepsize underflow in rkqs" << endl;*/ break;} 
        if(isnan(dt) || isinf(dt)) break;
    }
    if(!fsal) f(t + dt, x2, dx2);

    if(errmax > ERRCON) dt=SAFETY*dt*pow(errmax,PGROW);
    else dt *= 5.0;
    return dt;
}
