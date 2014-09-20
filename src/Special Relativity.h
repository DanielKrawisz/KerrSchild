#pragma once

//Note: None of these functions is actually used. 

//Boost a vector 'v' in the direction of a spacial vector 'x'.
// Set the boosted vector to 'newv'. 
void boost(double *newv, double *v, double *x);

//Set 'newm' as the boosted form of a matrix 'm' boosted in direction 'x'. 
void boostmtrx(double *newm, double *m, double *x);

//This returns the boost matrix but in the form of a symmetric tensor,
//  which means that only 10 numbers will be used rather than 16. 
void returntheboostmatrixasasymmetric2tensor(double *newm, double *v);

//Sets vdown as the dual vector to vup. 
void etadual(double vdown[], double vup[]);
