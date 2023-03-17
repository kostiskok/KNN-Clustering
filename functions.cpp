#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

#include "functions.h"
#include "data_handling.h"


using namespace std;

//extend mod for negative numbers
long mod(long D, long d){

    // if (d == 0){
    //     return D;
    // }

    if (D < 0){
        return (-D)%d;
    }
    else{
        return D%d;
    }

}

//L2 metric
double L2(pointType* p1, pointType* p2, int dims){
    double distance=0;
    for(int i=0;i<dims;i++){
        distance +=(p1->coords[i]-p2->coords[i])*(p1->coords[i]-p2->coords[i]);
    }
    return sqrt(distance);
}

double singlePointsL2(double p1, double p2){
    return sqrt((p1-p2)*(p1-p2));
}

//hamming distance between two numbers
int hammingDistance(int n1, int n2){
    int x = n1 ^ n2;
    int setBits = 0;
 
    while (x > 0) {
        setBits += x & 1;
        x >>= 1;
    }
 
    return setBits;
}

void add_x_values(vector<pointType>& curves){

    for (int c = 0; c < curves.size(); c++){

        int num_of_values = curves[c].coords.size();
        for (int i = num_of_values; i >= 1; i--){
            curves[c].coords.insert(curves[c].coords.begin()+i, i);
        }

    }

}

double min3(double p1, double p2, double p3){
    return min(min(p1,p2),p3);
}

int longestCurve(vector<pointType> curves){
    int max=0;
    for(int i=0;i<curves.size();i++){
        if(curves[i].coords.size()>max){
            max=curves[i].coords.size();
        }
    }
    //cout<<max<<endl;
    return max;
}

double maxValue(vector<pointType> curves){
    double max=numeric_limits<double>::min();
    for(int i=0;i<curves.size();i++){
        for(int j=0;j<curves[i].coords.size();j++){
            if(curves[i].coords[j]>max){
                max=curves[i].coords[j];
            }
        }
    }
    //cout<<max<<endl;
    return max;
}


void padding(vector<pointType>& curves){
    int paddingSize=longestCurve(curves);
    double paddingNumber=maxValue(curves);
    for(int i=0;i<curves.size();i++){
        for(int j=curves[i].coords.size();j<paddingSize;j++){
            curves[i].coords.push_back(paddingNumber);
        }
    }
}


double discreteFrechet(pointType curve1, pointType curve2){
    double answer;
    double **array=new double*[curve1.coords.size()];
    for(int i=0;i<curve1.coords.size();i++){
        array[i]=new double[curve2.coords.size()];
    }
    for(int i=1;i<curve1.coords.size();i+=2){
        for(int j=1;j<curve2.coords.size();j+=2){
            if(i==1 && j==1){
                //cout<<"1"<<endl;
                array[i][j]=L2(&curve1,&curve2,1);
            }
            else if(i==1){
                //cout<<"2"<<endl;
                array[i][j]=max(array[i][j-2],singlePointsL2(curve1.coords[i],curve2.coords[j])+singlePointsL2(curve1.coords[i-1],curve2.coords[j-1]));
            }
            else if(j==1){
                //cout<<"3"<<endl;
                array[i][j]=max(array[i-2][j],singlePointsL2(curve1.coords[i],curve2.coords[j])+singlePointsL2(curve1.coords[i-1],curve2.coords[j-1]));
            }
            else{
                //cout<<"4"<<endl;
                array[i][j]=max(min3(array[i-2][j],array[i][j-2],array[i-2][j-2]),singlePointsL2(curve1.coords[i],curve2.coords[j])+singlePointsL2(curve1.coords[i-1],curve2.coords[j-1]));
            }
        }
    }
    answer=array[curve1.coords.size()-1][curve2.coords.size()-1];

    for(int i=0;i<curve1.coords.size();i++){
        delete[] array[i];
    }
    delete[] array;
    //cout<<answer<<endl;
    return answer;
}

inline double snap(double coord,double delta,double t){
    //return floor((coord-t)/delta + 0.5)*delta + t;
    return floor(((coord + 0.5)/delta)*delta + t);
}
inline double snapContinuous(double coord,double delta,double t){
    return floor((coord+t)/delta)*delta;
}

pointType snapping2D(pointType curve, double delta, double* t){
    pointType snappedCurve;
    double x,y, prev_x, prev_y;

    //for (int i = 0; i < curve.coords.size()-1; i++){
        //curve.coords[i] += t[i];
    //}
    prev_x = snap(curve.coords[0],delta,t[0]);
    prev_y = snap(curve.coords[1],delta,t[1]);
    snappedCurve.coords.push_back(prev_x);
    snappedCurve.coords.push_back(prev_y);
    for(int i=2;i<curve.coords.size()-1;i+=2){
       x=snap(curve.coords[i],delta,t[i]);
       y=snap(curve.coords[i+1],delta,t[i+1]);
       if(x!=prev_x || y!=prev_y){
           snappedCurve.coords.push_back(x);
           snappedCurve.coords.push_back(y);
       }
       prev_x = x;
       prev_y = y;
    }
    return snappedCurve;
}

pointType snapping1D(pointType curve, double delta, double* t){
    pointType snappedCurve;
    double x,prev_x;
    //for (int i = 0; i < curve.coords.size()-1; i++){
        //curve.coords[i] += t[i];
    //}
    prev_x = snapContinuous(curve.coords[0],delta,t[0]);
    snappedCurve.coords.push_back(prev_x);
    for(int i=1;i<curve.coords.size();i++){
       x=snapContinuous(curve.coords[i],delta,t[i]);
       if(x!=prev_x){
           snappedCurve.coords.push_back(x);
       }
       prev_x = x;
    }
    return snappedCurve;



}

pointType filtering(pointType curve, double e=0.5){
    //for(int i=1;i<curve.coords.size()-1;i++){
    for(int i=curve.coords.size()-1;i > 0;i--){
        if(abs(curve.coords[i]-curve.coords[i-1])<e && abs(curve.coords[i]-curve.coords[i+1])<e ){
            curve.coords.erase(curve.coords.begin()+i);
        }
    }
    return curve;
}

double getMinDistance(vector<pointType> centroids,char metric){
    double min=numeric_limits<double>::max();
    double distance;
    for(int i=0;i<centroids.size();i++){
        for(int j=i+1;j<centroids.size();j++){
            if(metric='f'){
                distance=discreteFrechet(centroids[i],centroids[j]);
            }
            else{
                distance=L2(&centroids[i],&centroids[j],centroids[0].coords.size());
            }
            if(distance<min){
                min=distance;
            }
        }
    }
    return distance;
}