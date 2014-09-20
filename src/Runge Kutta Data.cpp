#include "Runge Kutta Data.h"
#include "Runge Kutta.h"

namespace RungeKutta{
    //Arrays describing various Runge-Kutta numerical methods. 
    double
        Euler_a[]            = {0.},
        Euler_b[]            = {1.},
        Euler_c[]            = {0.},
        Euler_d[]            = {0},
        Midpoint_a[]         = {2./3.},
        Midpoint_b[]         = {.25,          .75},
        Midpoint_c[]         = {2./3.},
        Midpoint_d[]         = {1.,0.},
        Kutta_a[]            = {.5,            0.,
                                -1.,           2.},
        Kutta_b[]            = {1/6.,          2/3.,          1/6.},
        Kutta_c[]            = {1/2.,          1.},
        Kutta_d[]            = {0.,            1.,            0.},
        BogackiSchampine_a[] = {1/2.,          0,             0,
                                0.,            3./4.,         0,
				2/9.,          1/3.,          4/9.},
        BogackiSchampine_b[] = {2/9.,          1/3.,          4/9.,          0},
        BogackiSchampine_c[] = {1/2.,          3./4.,         1.},
        BogackiSchampine_d[] = {7/24.,         1/4.,          1/3.,          1/8.},
        Fehlberg_a[]         = {1/4.,          0,             0,             0,             0,
                                3/32.,         9/32.,         0,             0,             0,
                                1932/2197.,    -7200./2197,   7296./2197,    0,             0,
                                439./216,      -8.,           3670/513.,     -845/4104.,    0,
                                -8/27.,        2.,            -3544/2565.,   1859/4104.,    -11/40},
        Fehlberg_d[]         = {25/216.,       0.,            1408/2565.,    2197/4104.,    -1/5.,         0},
        Fehlberg_c[]         = {1/4.,          3/8.,          12/13.,        1,             1/2.},
        Fehlberg_b[]         = {16/135.,       0.,            6656/12825.,   28561/56430.,  -9/50.,        2/55.},
        CashKarp_a[]         = {1/5.,          0.,            0.,            0.,            0.,
                                3/40.,         9/40.,         0.,            0.,            0.,
                                3/10.,         -9/10.,        6/5.,          0.,            0.,
                                -11/54.,       5/2.,          -70./27,       35/27.,        0.,
                                1631/55296.,   175/512.,      575/13824.,    44275/110592., 253/4096.},
        CashKarp_b[]         = {37/378.,       0.,            250/621.,      125./594,      0.,            512/1771.},
        CashKarp_c[]         = {1/5.,          3/10.,         3/5.,          1.,            7/8.},
        CashKarp_d[]         = {2825/27648.,   0.,            18575/48384.,  13525/55296.,  277/14336.,    1/4.},
        DormandPrince_a[]    = {1/5.,          0.,            0.,            0.,            0.,            0.,
                                3/40.,         9/40.,         0.,            0.,            0.,            0.,
                                44/45.,        -56/15.,       32/9.,         0.,            0.,            0.,
                                19372/6561.,   -25360/2187.,  64448/6561.,   -212/729.,     0.,            0.,
                                9017/3168.,    -355/33.,      46732/5247.,   49/176.,      -5103/18656.,   0.,
                                35/384.,       0.,            500/1113.,     125/192.,     -2187/6784.,    11./84},
        DormandPrince_b[]    = {35/384.,       0.,            500/1113.,     125/192.,     -2187/6784.,    11./84,        0},
        DormandPrince_c[]    = {1/5.,          3/10.,         4/5.,          8/9.,         1.,             1.},
        DormandPrince_d[]    = {5179/57600.,   0.,            7571/16695.,   393/640.,     -92097/339200., 187/2100.,     1/40.};

	//The arrays go into objects that are designed to implement the Runge-Kutta method. 
    RungeKuttaForm 
        *Euler			= new RungeKuttaForm(1, 1, 0, Euler_a, Euler_b, Euler_c, Euler_d),
        *Midpoint		= new RungeKuttaForm(2, 2, 0, Midpoint_a, Midpoint_b, Midpoint_c, Midpoint_d),
        *Kutta			= new RungeKuttaForm(3, 3, 0, Kutta_a, Kutta_b, Kutta_c, Kutta_d),
        *BogackiSchampine	= new RungeKuttaForm(4, 4, 1, BogackiSchampine_a, BogackiSchampine_b, BogackiSchampine_c, BogackiSchampine_d),
        *Fehlberg		= new RungeKuttaForm(6, 5, 0, Fehlberg_a, Fehlberg_b, Fehlberg_c, Fehlberg_d),
        //TODO: The Cash-Karp method just doesn't seem to work! The others do, however.
        //NOTE: I may have fixed the problem. There was a typo in the book I got it from.
        *CashKarp		= new RungeKuttaForm(6, 5, 0, CashKarp_a, CashKarp_b, CashKarp_c, CashKarp_d),
        *DormandPrince		= new RungeKuttaForm(7, 5, 1, DormandPrince_a, DormandPrince_b, DormandPrince_c, DormandPrince_d);
}
