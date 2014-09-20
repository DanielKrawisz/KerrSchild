#pragma once
#include <cmath>

//Just some simple math functions that I need. 

//The greater and lesser of two doubles. 
inline double dmax(double x, double y) {if (x>y) return x; else return y;}
inline double dmin(double x, double y) {if (x>y) return y; else return x;}

//trig functions. 
inline double cot(double x) {return 1/tan(x);}
inline double csc(double x) {return 1/sin(x);}
inline double sec(double x) {return 1/cos(x);}

inline int signof(double x) {return (x>0?1:-1);}
//inline bool isinf(double x) {return (x - 2 == x);}
