#ifndef LSH_FRECHET_H

#define LSH_FRECHET_H

#include <vector>

#include "HashTable_LSH.h"
#include "types.h"

class LSH_Frechet{

    private:
        std::vector<pointType>& curves;
        std::vector<pointType>* grid_curves;
        std::vector<pointType>* vector_x;

        HashTable_LSH** hashTables;

        int method; //0 for discrete, 1 for continuous

        const int k; //k h_i functions in the hash table
        const int L; //number of hash tables
        const double delta; //size of grid vertex
        const int w; //used by h_i (best values [4-6])

        double** t_value;

        static double continuous_frechet(pointType, pointType);

        void concat(pointType&, pointType);

    public:
        LSH_Frechet(std::vector<pointType>&, int, int, int, double, std::string);
        ~LSH_Frechet();

        void nearestNeighbours(pointType*, int, resultNN&);
        void rangeSearch(pointType*, int, resultR&);
};

#endif