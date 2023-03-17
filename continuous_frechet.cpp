#include "LSH_Frechet.h"
#include "fred_files/frechet.hpp"
#include "fred_files/curve.hpp"
#include "fred_files/point.hpp"

using namespace std;

double LSH_Frechet::continuous_frechet(pointType curve1, pointType curve2){
    
    //prepare the structures:
    //for the first curve
    Point** point1;
    point1 = new Point*[curve1.coords.size()];

    for(int i = 0; i < curve1.coords.size(); i++){
        point1[i] = new Point(1);
        (*point1[i]).set(0, curve1.coords[i]);
    }

    Points points1(curve1.coords.size(), 1);
    for (int i = 0; i < curve1.coords.size(); i++){
        points1.add(*point1[i]);
    }
    Curve c1(points1, curve1.id);

    //for the second curve
    Point** point2;
    point2 = new Point*[curve2.coords.size()];

    for(int i = 0; i < curve2.coords.size(); i++){
        point2[i] = new Point(1);
        (*point2[i]).set(0, curve2.coords[i]);
    }

    Points points2(curve2.coords.size(), 1);
    for (int i = 0; i < curve2.coords.size(); i++){
        points2.add(*point2[i]);
    }
    Curve c2(points2, curve2.id);

    //and calculate distance between them
    Frechet::Continuous::Distance dist = Frechet::Continuous::distance(c1, c2);

    //free memory
    for(int i = 0; i < curve1.coords.size(); i++){
        delete point1[i];
    }
    delete[] point1;

    for(int i = 0; i < curve2.coords.size(); i++){
        delete point2[i];
    }
    delete[] point2;
    //and return the distance
    return dist.value;
    
}