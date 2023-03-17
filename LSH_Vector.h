#ifndef LSH_VECTOR_H

#define LSH_VECTOR_H

#include <vector>

#include "HashTable_LSH.h"
#include "types.h"

class LSH_Vector{

    //pointer to function for the metric function
    typedef double (*metricFunction)(pointType*, pointType*, int);

    private:
        std::vector<pointType>& points; //reference to points, so as to
            //free all the memory already allocated before finishing

        HashTable_LSH** hashTables;

        const int k; //k h() per hash table
        const int L; //L hash tables
        const int d; //d dimensions per point/query
        const int w; //w "window"
        const metricFunction metric; //metric function

    public:
        LSH_Vector(std::vector<pointType>&, int, int, int, int, metricFunction);
        ~LSH_Vector();

        void nearestNeighbours(pointType*, int, resultNN&);
        void rangeSearch(pointType*, int, resultR&);

};

#endif