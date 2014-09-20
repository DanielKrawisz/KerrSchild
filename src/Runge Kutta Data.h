#pragma once

#include "Runge Kutta.h"

namespace RungeKutta{
    //These are arrays that contain the values of various Runge-Kutta 
    //  numerical methods. 
    extern double
        Euler_a[], Euler_b[], Euler_c[], Euler_d[],
        Midpoint_a[], Midpoint_b[], Midpoint_c[], Midpoint_d[],
        Kutta_a[], Kutta_b[], Kutta_c[], Kutta_d[],
        BogackiSchampine_a[], BogackiSchampine_b[], BogackiSchampine_c[], BogackiSchampine_d[], 
        Fehlberg_a[], Fehlberg_b[], Fehlberg_c[], Fehlberg_d[], 
        CashKarp_a[], CashKarp_b[], CashKarp_c[], CashKarp_d[], 
        DormandPrince_a[], DormandPrince_b[], DormandPrince_c[], DormandPrince_d[];
    //Objects that run different Runge Kutta methods. 
    extern RungeKuttaForm *Euler, *Midpoint, *Kutta, *BogackiSchampine, *Fehlberg, *CashKarp, *DormandPrince;
}
