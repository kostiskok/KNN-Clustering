#ifndef FUNCTIONS_H

#define FUNCTIONS_H

#include "types.h"

long mod(long, long);

double L2(pointType* ,pointType* , int);

double singlePointsL2(double , double );

int hammingDistance(int , int);

void add_x_values(std::vector<pointType>&);

double min3(double, double, double);

int longestCurve(std::vector<pointType>);

double maxValue(std::vector<pointType>);

void padding(std::vector<pointType>&);

double discreteFrechet(pointType, pointType);

pointType snapping2D(pointType, double, double*);

pointType snapping1D(pointType , double, double*);

pointType filtering(pointType, double);

double getMinDistance(std::vector<pointType>,char);

#endif