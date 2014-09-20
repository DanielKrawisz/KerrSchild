//TODO: the metric functions ought to return the metric with all 16 numbers. 
//TODO: Check what coordinate system we should be in when the file is loaded.

#pragma once

#define PIXELNOTIFICATION		1000
#define STEPNOTIFICATION		10000
#define DISTANCERANGE			10.0
#define MAXOUTERADJUSTMENTS		20
#define INPUT				"ksTest.txt"
#define MAXSTEPS			100000

#define ITEM                    15
#define IND                     "   "

void RunUniverse(const char *file);
