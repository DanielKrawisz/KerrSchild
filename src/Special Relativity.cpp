#include "Special Relativity.h"
#include "Extra Math.h"

//This also needs to be tested. 
void boost(double *newv, double *v, double *x)
{
    static double b, g, gg;
    b = v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
    if(isinf(1/b) || isnan(1/b)) 
    {
        for(int i = 0; i < 4; i++) newv[i] = x[i];
    }
    else
    {
        g = 1/sqrt(1 - b);
        gg = (g - 1)/b;

        newv[0] = g*x[0] - g*v[1]*x[1] - g*v[2]*x[2] - g*v[3]*x[3];
        newv[1] = -g*v[1]*x[0] + (1 + gg*v[1]*v[1])*x[1] + gg*v[1]*v[2]*x[2] + gg*v[1]*v[3]*x[3];
        newv[2] = -g*v[2]*x[0] + gg*v[1]*v[2]*x[1] + (1 + gg*v[2]*v[2])*x[2] + gg*v[2]*v[3]*x[3];
        newv[3] = -g*v[3]*x[0] + gg*v[1]*v[3]*x[1] + gg*v[2]*v[3]*x[2] + (1 + gg*v[3]*v[3])*x[3];
    }
}

//Set 'newm' as the boosted form of a matrix 'm' boosted in direction 'x'.
void boostmtrx(double *newm, double *m, double *x)
{
    static double b, g, gg;
    b = x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
    if(isinf(1/b) || isnan(1/b)) 
    {
        for(int i = 0; i < 16; i++) newm[i] = m[i];
    }
    else
    {
        g = 1/sqrt(1 - b);
        gg = (g - 1)/b;

        newm[0] = g*(m[0] + m[5]*x[1]*x[1] - m[2]*x[2] - m[8]*x[2] + m[10]*x[2]*x[2] - m[3]*x[3] - m[12]*x[3] + m[11]*x[2]*x[3] + m[14]*x[2]*x[3] + m[15]*x[3]*x[3] + x[1]*(-m[1] - m[4] + m[6]*x[2] + m[9]*x[2] + m[7]*x[3] + m[13]*x[3]));
        newm[1] = g*(g*x[1]*(-m[0] + m[4]*x[1] + m[8]*x[2] + m[12]*x[3]) - (1 + gg*x[1]*x[1])*(-m[1] + m[5]*x[1] + m[9]*x[2] + m[13]*x[3]) - gg*x[1]*x[2]*(-m[2] + m[6]*x[1] + m[10]*x[2] + m[14]*x[3]) - gg*x[1]*x[3]*(-m[3] + m[7]*x[1] + m[11]*x[2] + m[15]*x[3]));
        newm[2] = g*(g*x[2]*(-m[0] + m[4]*x[1] + m[8]*x[2] + m[12]*x[3]) - gg*x[1]*x[2]*(-m[1] + m[5]*x[1] + m[9]*x[2] + m[13]*x[3]) - (1 + gg*x[2]*x[2])*(-m[2] + m[6]*x[1] + m[10]*x[2] + m[14]*x[3]) - gg*x[2]*x[3]*(-m[3] + m[7]*x[1] + m[11]*x[2] + m[15]*x[3]));
        newm[3] = g*(g*x[3]*(-m[0] + m[4]*x[1] + m[8]*x[2] + m[12]*x[3]) - gg*x[1]*x[3]*(-m[1] + m[5]*x[1] + m[9]*x[2] + m[13]*x[3]) - gg*x[2]*x[3]*(-m[2] + m[6]*x[1] + m[10]*x[2] + m[14]*x[3]) - (-m[3] + m[7]*x[1] + m[11]*x[2] + m[15]*x[3])*(1 + gg*x[3]*x[3]));
        newm[4] = g*(-(g*m[0]*x[1]) + m[4]*(1 + gg*x[1]*x[1]) + gg*m[8]*x[1]*x[2] + gg*m[12]*x[1]*x[3] - x[1]*(m[5] - g*m[1]*x[1] + gg*m[5]*x[1]*x[1] + gg*m[9]*x[1]*x[2] + gg*m[13]*x[1]*x[3]) - x[2]*(m[6] - g*m[2]*x[1] + gg*m[6]*x[1]*x[1] + gg*m[10]*x[1]*x[2] + gg*m[14]*x[1]*x[3]) - x[3]*(m[7] - g*m[3]*x[1] + gg*m[7]*x[1]*x[1] + gg*m[11]*x[1]*x[2] + gg*m[15]*x[1]*x[3]));
        newm[5] = -(g*x[1]*(m[4] - g*m[0]*x[1] + gg*m[4]*x[1]*x[1] + gg*m[8]*x[1]*x[2] + gg*m[12]*x[1]*x[3])) + (1 + gg*x[1]*x[1])*(m[5] - g*m[1]*x[1] + gg*m[5]*x[1]*x[1] + gg*m[9]*x[1]*x[2] + gg*m[13]*x[1]*x[3]) + gg*x[1]*x[2]*(m[6] - g*m[2]*x[1] + gg*m[6]*x[1]*x[1] + gg*m[10]*x[1]*x[2] + gg*m[14]*x[1]*x[3]) + gg*x[1]*x[3]*(m[7] - g*m[3]*x[1] + gg*m[7]*x[1]*x[1] + gg*m[11]*x[1]*x[2] + gg*m[15]*x[1]*x[3]);
        newm[6] = -(g*x[2]*(m[4] - g*m[0]*x[1] + gg*m[4]*x[1]*x[1] + gg*m[8]*x[1]*x[2] + gg*m[12]*x[1]*x[3])) + gg*x[1]*x[2]*(m[5] - g*m[1]*x[1] + gg*m[5]*x[1]*x[1] + gg*m[9]*x[1]*x[2] + gg*m[13]*x[1]*x[3]) + (1 + gg*x[2]*x[2])*(m[6] - g*m[2]*x[1] + gg*m[6]*x[1]*x[1] + gg*m[10]*x[1]*x[2] + gg*m[14]*x[1]*x[3]) + gg*x[2]*x[3]*(m[7] - g*m[3]*x[1] + gg*m[7]*x[1]*x[1] + gg*m[11]*x[1]*x[2] + gg*m[15]*x[1]*x[3]);
        newm[7] = -(g*x[3]*(m[4] - g*m[0]*x[1] + gg*m[4]*x[1]*x[1] + gg*m[8]*x[1]*x[2] + gg*m[12]*x[1]*x[3])) + gg*x[1]*x[3]*(m[5] - g*m[1]*x[1] + gg*m[5]*x[1]*x[1] + gg*m[9]*x[1]*x[2] + gg*m[13]*x[1]*x[3]) + gg*x[2]*x[3]*(m[6] - g*m[2]*x[1] + gg*m[6]*x[1]*x[1] + gg*m[10]*x[1]*x[2] + gg*m[14]*x[1]*x[3]) + (m[7] - g*m[3]*x[1] + gg*m[7]*x[1]*x[1] + gg*m[11]*x[1]*x[2] + gg*m[15]*x[1]*x[3])*(1 + gg*x[3]*x[3]);
        newm[8] = g*(m[8] - m[9]*x[1] - m[10]*x[2] - m[11]*x[3] + g*x[2]*(-m[0] + m[1]*x[1] + m[2]*x[2] + m[3]*x[3]) - gg*x[2]*(m[5]*x[1]*x[1] - m[8]*x[2] + m[10]*x[2]*x[2] - m[12]*x[3] + m[11]*x[2]*x[3] + m[14]*x[2]*x[3] + m[15]*x[3]*x[3] + x[1]*(-m[4] + (m[6] + m[9])*x[2] + (m[7] + m[13])*x[3])));
        newm[9] = -(g*x[1]*(m[8] - g*m[0]*x[2] + gg*x[2]*(m[4]*x[1] + m[8]*x[2] + m[12]*x[3]))) + (1 + gg*x[1]*x[1])*(m[9] - g*m[1]*x[2] + gg*x[2]*(m[5]*x[1] + m[9]*x[2] + m[13]*x[3])) + gg*x[1]*x[2]*(m[10] - g*m[2]*x[2] + gg*x[2]*(m[6]*x[1] + m[10]*x[2] + m[14]*x[3])) + gg*x[1]*x[3]*(m[11] - g*m[3]*x[2] + gg*x[2]*(m[7]*x[1] + m[11]*x[2] + m[15]*x[3]));
        newm[10] = -(g*x[2]*(m[8] - g*m[0]*x[2] + gg*x[2]*(m[4]*x[1] + m[8]*x[2] + m[12]*x[3]))) + gg*x[1]*x[2]*(m[9] - g*m[1]*x[2] + gg*x[2]*(m[5]*x[1] + m[9]*x[2] + m[13]*x[3])) + (1 + gg*x[2]*x[2])*(m[10] - g*m[2]*x[2] + gg*x[2]*(m[6]*x[1] + m[10]*x[2] + m[14]*x[3])) + gg*x[2]*x[3]*(m[11] - g*m[3]*x[2] + gg*x[2]*(m[7]*x[1] + m[11]*x[2] + m[15]*x[3]));
        newm[11] = -(g*x[3]*(m[8] - g*m[0]*x[2] + gg*x[2]*(m[4]*x[1] + m[8]*x[2] + m[12]*x[3]))) + gg*x[1]*x[3]*(m[9] - g*m[1]*x[2] + gg*x[2]*(m[5]*x[1] + m[9]*x[2] + m[13]*x[3])) + gg*x[2]*x[3]*(m[10] - g*m[2]*x[2] + gg*x[2]*(m[6]*x[1] + m[10]*x[2] + m[14]*x[3])) + (1 + gg*x[3]*x[3])*(m[11] - g*m[3]*x[2] + gg*x[2]*(m[7]*x[1] + m[11]*x[2] + m[15]*x[3]));
        newm[12] = g*(m[12] - m[13]*x[1] - m[14]*x[2] - m[15]*x[3] + g*x[3]*(-m[0] + m[1]*x[1] + m[2]*x[2] + m[3]*x[3]) - gg*x[3]*(m[5]*x[1]*x[1] - m[8]*x[2] + m[10]*x[2]*x[2] - m[12]*x[3] + m[11]*x[2]*x[3] + m[14]*x[2]*x[3] + m[15]*x[3]*x[3] + x[1]*(-m[4] + (m[6] + m[9])*x[2] + (m[7] + m[13])*x[3])));
        newm[13] = -(g*x[1]*(m[12] - g*m[0]*x[3] + gg*x[3]*(m[4]*x[1] + m[8]*x[2] + m[12]*x[3]))) + (1 + gg*x[1]*x[1])*(m[13] - g*m[1]*x[3] + gg*x[3]*(m[5]*x[1] + m[9]*x[2] + m[13]*x[3])) + gg*x[1]*x[2]*(m[14] - g*m[2]*x[3] + gg*x[3]*(m[6]*x[1] + m[10]*x[2] + m[14]*x[3])) + gg*x[1]*x[3]*(m[15] - g*m[3]*x[3] + gg*x[3]*(m[7]*x[1] + m[11]*x[2] + m[15]*x[3]));
        newm[14] = -(g*x[2]*(m[12] - g*m[0]*x[3] + gg*x[3]*(m[4]*x[1] + m[8]*x[2] + m[12]*x[3]))) + gg*x[1]*x[2]*(m[13] - g*m[1]*x[3] + gg*x[3]*(m[5]*x[1] + m[9]*x[2] + m[13]*x[3])) + (1 + gg*x[2]*x[2])*(m[14] - g*m[2]*x[3] + gg*x[3]*(m[6]*x[1] + m[10]*x[2] + m[14]*x[3])) + gg*x[2]*x[3]*(m[15] - g*m[3]*x[3] + gg*x[3]*(m[7]*x[1] + m[11]*x[2] + m[15]*x[3]));
        newm[15] = -(g*x[3]*(m[12] - g*m[0]*x[3] + gg*x[3]*(m[4]*x[1] + m[8]*x[2] + m[12]*x[3]))) + gg*x[1]*x[3]*(m[13] - g*m[1]*x[3] + gg*x[3]*(m[5]*x[1] + m[9]*x[2] + m[13]*x[3])) + gg*x[2]*x[3]*(m[14] - g*m[2]*x[3] + gg*x[3]*(m[6]*x[1] + m[10]*x[2] + m[14]*x[3])) + (1 + gg*x[3]*x[3])*(m[15] - g*m[3]*x[3] + gg*x[3]*(m[7]*x[1] + m[11]*x[2] + m[15]*x[3]));
    }
}

//This returns the boost matrix but in the form of a symmetric tensor,
//  which means that only 10 numbers will be used rather than 16. 
void returntheboostmatrixasasymmetric2tensor(double *m, double *v)
{
    static double b, g, gg;
    b = v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
    if(isinf(1/b) || isnan(1/b)) 
    {
        m[1] = 0; m[2] = 0; m[3] = 0; m[6] = 0; m[5] = 0; m[8] = 0;
        m[0] = 1; m[4] = 1; m[7] = 1; m[9] = 1;
    }
    else
    {
        g = 1/sqrt(1 - b);
        gg = (g - 1)/b;

        m[0] = g;
        m[1] = -v[1] * g;
        m[2] = -v[2] * g;
        m[3] = -v[3] * g;
        m[4] = 1 + gg * v[1] * v[1];
        m[5] = gg * v[1] * v[2];
        m[6] = gg * v[1] * v[3];
        m[7] = 1 + gg * v[2] * v[2];
        m[8] = gg * v[2] * v[3];
        m[9] = 1 + gg * v[3] * v[3];
    }
}

//Sets vdown as the dual vector to vup. 
void etadual(double vdown[], double vup[])
{
    vdown[0] = vup[0]; vdown[1] = -vup[1]; vdown[2] = -vup[2]; vdown[3] = -vup[3];
}
