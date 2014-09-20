#pragma once

//All these functions assume four dimensions since that is the number
//  of dimensions general relativity takes place in. 

//Take two 4-vectors and a 4x4 metric and return their inner product. 
double InnerProduct(double* v1, double* v2, double* metric);

//Set a vector to be the product of a 4x4 matrix and a 4-vector. 
void MatrixDotVector(double* newvec, double* mtrx, double* vec);

//Transpose a 4x4 matrix
void Transpose(double *mtrx);

//set 4-vector 'range' to be a copy of 4-vector 'domain'. 
void CopyVector(double* range, double* domain);

//Set 'inverse' to be the inverse of 4x4 matrix 'm'
void InverseMatrix(double* inverse, double* m);

//Given a set of four 4-vectors and a metric, 
//  turn them into an orthonormal set. 
void Orthonormalize(double* mtrx, double* metric);
