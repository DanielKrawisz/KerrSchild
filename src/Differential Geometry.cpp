#include "Differential Geometry.h"
#include <cmath>

using namespace std;

//Take two 4-vectors and a 4x4 metric and return their inner product. 
double InnerProduct(double* v1, double* v2, double* metric)
{
    static double result;
    static double downvec[4];
    static int i;
    result = 0.0;

    MatrixDotVector(downvec, metric, v2);
    for(i=0;i<4;i++) result += v1[i]*downvec[i];

    return result;
}

//Set a vector to be the product of a 4x4 matrix and a 4-vector. 
void MatrixDotVector(double* newvec, double* mtrx, double* vec)
{
    static double temp[4];
    temp[0] =  mtrx[0]*vec[0] +  mtrx[1]*vec[1] +  mtrx[2]*vec[2] +  mtrx[3]*vec[3];
    temp[1] =  mtrx[4]*vec[0] +  mtrx[5]*vec[1] +  mtrx[6]*vec[2] +  mtrx[7]*vec[3];
    temp[2] =  mtrx[8]*vec[0] +  mtrx[9]*vec[1] + mtrx[10]*vec[2] + mtrx[11]*vec[3];
    temp[3] = mtrx[12]*vec[0] + mtrx[13]*vec[1] + mtrx[14]*vec[2] + mtrx[15]*vec[3];
    newvec[0] = temp[0]; newvec[1] = temp[1]; newvec[2] = temp[2]; newvec[3] = temp[3];
}

//Transpose a 4x4 matrix
void Transpose(double *mtrx)
{
    int i, j;
    double tmp;
    for(i = 0; i < 4; i++)
        for(j = i + 1; j < 4; j++)
        {
            tmp = mtrx[4*i+j];
            mtrx[4*i+j] = mtrx[4*j+i];
            mtrx[4*j+i] = tmp;
        }
}

//set 4-vector 'range' to be a copy of 4-vector 'domain'.
void CopyVector(double* range, double* domain)
{
    for(int i = 0; i < 4; i++) range[i] = domain[i];
}

//Set 'inverse' to be the inverse of 4x4 matrix 'm'
void InverseMatrix(double* inverse, double* m)
{
    static double det;
    det = m[1]*m[11]*m[14]*m[4] - m[1]*m[10]*m[15]*m[4] - m[11]*m[13]*m[2]*m[4] + m[10]*m[13]*m[3]*m[4] - m[0]*m[11]*m[14]*m[5] + m[0]*m[10]*m[15]*m[5] + m[11]*m[12]*m[2]*m[5] - m[10]*m[12]*m[3]*m[5] - m[1]*m[11]*m[12]*m[6] + m[0]*m[11]*m[13]*m[6] + m[1]*m[10]*m[12]*m[7] - m[0]*m[10]*m[13]*m[7] - m[15]*m[2]*m[5]*m[8] + m[14]*m[3]*m[5]*m[8] + m[1]*m[15]*m[6]*m[8] - m[13]*m[3]*m[6]*m[8] - m[1]*m[14]*m[7]*m[8] + m[13]*m[2]*m[7]*m[8] + m[15]*m[2]*m[4]*m[9] - m[14]*m[3]*m[4]*m[9] - m[0]*m[15]*m[6]*m[9] + m[12]*m[3]*m[6]*m[9] + m[0]*m[14]*m[7]*m[9] - m[12]*m[2]*m[7]*m[9];

    inverse[0] = (-(m[11]*m[14]*m[5]) + m[10]*m[15]*m[5] + m[11]*m[13]*m[6] - m[10]*m[13]*m[7] - m[15]*m[6]*m[9] + m[14]*m[7]*m[9])/det;
    inverse[1] = (m[1]*m[11]*m[14] - m[1]*m[10]*m[15] - m[11]*m[13]*m[2] + m[10]*m[13]*m[3] + m[15]*m[2]*m[9] - m[14]*m[3]*m[9])/det;
    inverse[2] = (-(m[15]*m[2]*m[5]) + m[14]*m[3]*m[5] + m[1]*m[15]*m[6] - m[13]*m[3]*m[6] - m[1]*m[14]*m[7] + m[13]*m[2]*m[7])/det;
    inverse[3] = (m[11]*m[2]*m[5] - m[10]*m[3]*m[5] - m[1]*m[11]*m[6] + m[1]*m[10]*m[7] + m[3]*m[6]*m[9] - m[2]*m[7]*m[9])/det;
    inverse[4] = (m[11]*m[14]*m[4] - m[10]*m[15]*m[4] - m[11]*m[12]*m[6] + m[10]*m[12]*m[7] + m[15]*m[6]*m[8] - m[14]*m[7]*m[8])/det;
    inverse[5] = (-(m[0]*m[11]*m[14]) + m[0]*m[10]*m[15] + m[11]*m[12]*m[2] - m[10]*m[12]*m[3] - m[15]*m[2]*m[8] + m[14]*m[3]*m[8])/det;
    inverse[6] = (m[15]*m[2]*m[4] - m[14]*m[3]*m[4] - m[0]*m[15]*m[6] + m[12]*m[3]*m[6] + m[0]*m[14]*m[7] - m[12]*m[2]*m[7])/det;
    inverse[7] = (-(m[11]*m[2]*m[4]) + m[10]*m[3]*m[4] + m[0]*m[11]*m[6] - m[0]*m[10]*m[7] - m[3]*m[6]*m[8] + m[2]*m[7]*m[8])/det;
    inverse[8] = (-(m[11]*m[13]*m[4]) + m[11]*m[12]*m[5] - m[15]*m[5]*m[8] + m[13]*m[7]*m[8] + m[15]*m[4]*m[9] - m[12]*m[7]*m[9])/det;
    inverse[9] = (-(m[1]*m[11]*m[12]) + m[0]*m[11]*m[13] + m[1]*m[15]*m[8] - m[13]*m[3]*m[8] - m[0]*m[15]*m[9] + m[12]*m[3]*m[9])/det;
    inverse[10] = (-(m[1]*m[15]*m[4]) + m[13]*m[3]*m[4] + m[0]*m[15]*m[5] - m[12]*m[3]*m[5] + m[1]*m[12]*m[7] - m[0]*m[13]*m[7])/det;
    inverse[11] = (m[1]*m[11]*m[4] - m[0]*m[11]*m[5] + m[3]*m[5]*m[8] - m[1]*m[7]*m[8] - m[3]*m[4]*m[9] + m[0]*m[7]*m[9])/det;
    inverse[12] = (m[10]*m[13]*m[4] - m[10]*m[12]*m[5] + m[14]*m[5]*m[8] - m[13]*m[6]*m[8] - m[14]*m[4]*m[9] + m[12]*m[6]*m[9])/det;
    inverse[13] = (m[1]*m[10]*m[12] - m[0]*m[10]*m[13] - m[1]*m[14]*m[8] + m[13]*m[2]*m[8] + m[0]*m[14]*m[9] - m[12]*m[2]*m[9])/det;
    inverse[14] = (m[1]*m[14]*m[4] - m[13]*m[2]*m[4] - m[0]*m[14]*m[5] + m[12]*m[2]*m[5] - m[1]*m[12]*m[6] + m[0]*m[13]*m[6])/det;
    inverse[15] = (-(m[1]*m[10]*m[4]) + m[0]*m[10]*m[5] - m[2]*m[5]*m[8] + m[1]*m[6]*m[8] + m[2]*m[4]*m[9] - m[0]*m[6]*m[9])/det;
}

//Given a set of four 4-vectors and a metric, 
//  turn them into an orthonormal set. 
void Orthonormalize(double* mtrx, double* metric)
{
    double ip;
    int i, j, k;
    double *v1, *v2;
    //This loop simply uses the Gram-Schmidt process.
    for(i = 0; i < 16; i+=4)
    {
        v1 = &(mtrx[i]); 
        for(j = 0; j < i; j+=4)
        {
            v2 = &(mtrx[j]); 
            ip = InnerProduct(v1, v2, metric);
            for(k = 0; k < 4; k++) mtrx[i+k] -= ip*mtrx[j+k];
        }

        ip = InnerProduct(v1, v1, metric);
        ip = pow(abs(ip), -.5);
        for(k = 0; k < 4; k++) mtrx[i+k] *= ip;
    }
}
