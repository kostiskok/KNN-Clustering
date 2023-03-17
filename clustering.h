#ifndef CLUSTERING_H

#define CLUSTERING_H

#include <vector>
#include "types.h"
#include "LSH_Vector.h"
#include "Hypercube.h"
#include "LSH_Frechet.h"
#include "data_handling.h"

class clustering{

    //pointer to function for the metric function
    typedef double (*metricFunction)(pointType*, pointType*, int);

    private:
        std::vector<pointType>& points; //reference to points, so as to
            //free all the memory already allocated before finishing

        std::vector<pointType> centroids;
        std::vector<std::vector<pointType*>> points_per_centroid;
        std::vector<double> silhouette;

        pointType temp_cluster;
        std::map<std::string, bool> isAssigned; //used by assign_LSH and assign_Hypercube

        char method; // c -> classic (Lloyd's), l -> LSH, h -> Hypercube, f-> Frechet
        char update; // v -> mean vector, f -> mean frechet

        double limit; //limit for update functions

        //parameters:
        const int K; // number of clusters/K-medians
        int d; // number of dimensions per point

        int L; // number of hashtables (for LSH)

        int M; // max number of elements to be checked (for Hypercube)
        int probes; // max number of probes to be checked (for Hypercube)

        int k; // k number of h() function / f() (for LSH and Hypercube)

        metricFunction metric;

        //different algorithms:
        LSH_Vector *lsh_v;
        Hypercube *hc;
        LSH_Frechet *lsh_f;

        std::vector<pointType> clusters;

        void initialize_common();
        void assignPoint(pointType&, pointType);
        double minimum_distance(pointType, int&, int);

        void k_means();

        void assignment_lloyds();
        void assignment_rangeSearch_LSH_Vector();
        void assignment_rangeSearch_Hypercube();
        void assignment_rangeSearch_LSH_Frechet();

        bool update_Vector();
        bool update_Frechet();

    public:
        clustering(std::vector<pointType>&, int, int, char);
        clustering(std::vector<pointType>&, int, int, int, int);
        clustering(std::vector<pointType>&, int, int, int, int, int);
        clustering(std::vector<pointType>&, int, int, int, int, double);
        ~clustering();

        void cluster(resultCL&, bool);
        
};

#endif