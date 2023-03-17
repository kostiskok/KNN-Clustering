#include <iostream>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <random>

#include "LSH_Frechet.h"
#include "functions.h"
#include "data_handling.h"

using namespace std;

LSH_Frechet::LSH_Frechet(vector<pointType>& c, int k, int l, int w, double d, string m): k(k), L(l), delta(d), w(w), curves(c){

    vector<pointType> curves_filt(c.size());

    if (m == "discrete"){
        method = 0;
    }
    else{
        method = 1;
    }

    //for every grid (aka L times), map each curve to every grid (aka snap them)
    //-> delta is common among all grids
    //-> t is random -> different grids

    //each curve originally has the same size
    int size = curves[0].coords.size();

    random_device r;
    default_random_engine gen(r());
    uniform_real_distribution<double> uni_d(0.0, delta);

    //calculate for each grid the random t vector (R^delta)
    t_value = new double*[L];
    for (int i = 0; i < L; i++){
        t_value[i] = new double[size];
        for (int j = 0; j < size; j++){
            t_value[i][j] = uni_d(gen);
        }
    }

    if (method == 1){
        
        for (int i = 0; i < curves.size(); i++){
            curves_filt[i] = filtering(curves[i], 1/2);
        }

    }

    //each curve is mapped to L curves, with different t
    grid_curves = new vector<pointType>[L];
    for (int i = 0; i < L; i++){

        if (method == 0){
            grid_curves[i].resize(curves.size());
        }
        else{
            grid_curves[i].resize(curves_filt.size());
        }

        for (int j = 0; j < curves.size(); j++){

            if (method == 0){
                grid_curves[i][j] = snapping2D(curves[j], delta, t_value[i]);
            }
            else{
                grid_curves[i][j] = snapping1D(curves_filt[j], delta, t_value[i]);
            }

        }

    }

    for (int i = 0; i < L; i++){
        padding(grid_curves[i]);
    }

    //then represent each grid curve to a vector (aka concat) 
    //this is the real vector x -- only discrete

    vector_x = new vector<pointType>[L];
    for (int i = 0; i < L; i++){
        vector_x[i].resize(curves.size());

        for (int j = 0; j < curves.size(); j++){

            vector_x[i][j].coords.resize(grid_curves[i][j].coords.size()/2);

            //concat(vector_x[i][j], grid_curves[i][j]);
            if (method == 0){
                for (int k = 0; k < grid_curves[i][j].coords.size(); k+=2){
                    vector_x[i][j].coords[k/2] = grid_curves[i][j].coords[k] + grid_curves[i][j].coords[k+1];
                }
            }
            else{
                for (int k = 0; k < grid_curves[i][j].coords.size(); k++){
                    vector_x[i][j].coords[k] + grid_curves[i][j].coords[k];
                }
            }

            //associate x with pointer to the grid_curve and to curve
            vector_x[i][j].orig_curve = &curves[j];
            vector_x[i][j].grid_curve = &grid_curves[i][j];
        }
    }

    //and now for each of the L hash tables, insert each vector x

    int tableSize = curves.size()/32+1; //Hash table's size: n/8 or n/16 -> n/32
    
    //Initialize L hash tables, and insert each vector_x in the right one
    //(aka vector_x[i] -> hashTable[i])
    hashTables = new HashTable_LSH*[L];
    
    for (int i = 0; i < L; i++){
        int coords_size = vector_x[i][0].coords.size();
        hashTables[i] = new HashTable_LSH(tableSize, k, coords_size, w);

        //!! insert not the curve but the real vector x
        for (int j = 0; j < curves.size(); j++){

            vector_x[i][j].hash_id = new long[L];

            hashTables[i]->insertPoint(&vector_x[i][j]);
        }

    }

}

//Free all memory that was previously allocated
LSH_Frechet::~LSH_Frechet(){

    for(int i = 0; i < L; i++){
        delete[] t_value[i];
    }
    delete[] t_value;

    delete[] grid_curves;

    for (int i = 0; i < L; i++){
        delete hashTables[i];
    }
    delete[] hashTables;

    for (int i = 0; i < L; i++){
        for (int j = 0; j < curves.size(); j++){
            delete[] vector_x[i][j].hash_id;
        }
    }
    delete[] vector_x;

}

//void LSH::nearestNeighbours(pointType query, resultNN& rNN){
void LSH_Frechet::nearestNeighbours(pointType* query, int N, resultNN& rNN){

    double dist, distTrue;
    int query_id, numOfHashTables = L;
    bool flag = false;

    vector<pointType*> bucket;
    vector<pair<double, string>> dTrue(N); //True best distance
    vector<pair<double, string>> bestDist(N); //for each q {best distance, best candidate}
    
    vector<pointType> grid_curve_query(numOfHashTables);
    vector<pointType> vector_x_query(numOfHashTables);

    //Initialiazation
    for (int i = 0; i < N; i++){
        bestDist[i] = {numeric_limits<double>::max(), "EMPTY"};
        dTrue[i] = {numeric_limits<double>::max(), "EMPTY"};
    }

    auto startTrue = chrono::high_resolution_clock::now();

    //Calculate true best distances (brute force) and store the distances (Not the points)
    for (int p = 0; p < curves.size(); p++){

        if (method == 0){
            distTrue = discreteFrechet(curves[p], *query);
        }
        else{
            distTrue = continuous_frechet(curves[p], *query);
        }
        dTrue.insert(upper_bound(dTrue.cbegin(),dTrue.cend(), make_pair(distTrue, curves[p].id)),make_pair(distTrue, curves[p].id));
        dTrue.pop_back();
    }

    auto stopTrue=chrono::high_resolution_clock::now();
    auto durationTrue=chrono::duration<double>(stopTrue-startTrue);

    auto start = chrono::high_resolution_clock::now();

    for (int i = 0; i < numOfHashTables; i++){

        if (method == 0){
            grid_curve_query[i] = snapping2D(*query, delta, t_value[i]);
        }
        else{
            pointType query_filt = filtering(*query, 1/2);
            grid_curve_query[i] = snapping1D(query_filt, delta, t_value[i]);
        }

        //concat(vector_x_query[i], grid_curve_query[i]);
        if (method == 0){
            vector_x_query[i].coords.resize(grid_curve_query[i].coords.size()/2);
            for (int k = 0; k < grid_curve_query[i].coords.size(); k += 2){
                vector_x_query[i].coords[k/2] = grid_curve_query[i].coords[k] + grid_curve_query[i].coords[k+1];
            }
        }
        else{
            vector_x_query[i].coords.resize(grid_curve_query[i].coords.size());
            for (int k = 0; k < grid_curve_query[i].coords.size(); k ++){
                vector_x_query[i].coords[k] = grid_curve_query[i].coords[k];
            }
        }

        vector_x_query[i].hash_id = new long[L];
        hashTables[i]->getBucket(&vector_x_query[i], bucket);
        query_id = vector_x_query[i].hash_id[i];

        flag = false;

        for (int j = 0; j < bucket.size(); j++){

            //Initially calculate distances from points in the bucket with same id from q (*)
            if (bucket[j]->hash_id[i] == query_id){

                flag = true;

                if (method == 0){
                    dist = discreteFrechet(*bucket[j]->orig_curve, *query);
                }
                else{
                    dist = continuous_frechet(*bucket[j]->orig_curve, *query);
                }

                for (int index = 0; index < N; index++){

                    //If point already in the bucket, ignore it
                    if (bucket[j]->orig_curve->id == bestDist[index].second){
                        break;
                    }

                    //insert sorted in bucket {distance from point_j, name of point_j}
                    if (dist < bestDist[index].first){
                        bestDist.insert(upper_bound(bestDist.begin(), bestDist.end(), make_pair(dist, bucket[j]->orig_curve->id)), make_pair(dist, bucket[j]->orig_curve->id));
                        bestDist.pop_back();
                        break;
                    }

                }

            }

        }

        //If no such points where found (*), then calculate the distance of each point in the bucket
        if (flag == false){

            for (int j = 0; j < bucket.size(); j++){

                if (method == 0){
                    dist = discreteFrechet(*bucket[j]->orig_curve, *query);
                }
                else{
                    dist = continuous_frechet(*bucket[j]->orig_curve, *query);
                }

                for (int index = 0; index < N; index++){

                    //If point already in the bucket, ignore it
                    if (bucket[j]->orig_curve->id == bestDist[index].second){
                        break;
                    }

                    //insert sorted in bucket {distance from point_j, name of point_j}
                    if (dist < bestDist[index].first){
                        bestDist.insert(upper_bound(bestDist.begin(), bestDist.end(), make_pair(dist, bucket[j]->orig_curve->id)), make_pair(dist, bucket[j]->orig_curve->id));
                        bestDist.pop_back();
                        break;
                    }

                }

            }

        }

    }

    for (int i = 0; i < numOfHashTables; i++){
        delete[] vector_x_query[i].hash_id;
    }

    auto stop=chrono::high_resolution_clock::now();
    auto duration=chrono::duration<double>(stop-start);

    rNN.Dist = bestDist;
    rNN.trueDist = dTrue;
    rNN.tDist = duration.count();
    rNN.tTrue = durationTrue.count();

}

void LSH_Frechet::rangeSearch(pointType* query, int R, resultR& rR){

    double dist;
    int query_id, numOfHashTables = L;

    vector<pointType*> bucket;
    vector<pointType*> dR; //points within r distance from each q

    vector<pointType> grid_curve_query(numOfHashTables);
    vector<pointType> vector_x_query(numOfHashTables);

    for (int i = 0; i < numOfHashTables; i++){

        if (method == 0){
            grid_curve_query[i] = snapping2D(*query, delta, t_value[i]);
        }
        else{
            pointType query_filt = filtering(*query, 1/2);
            grid_curve_query[i] = snapping1D(*query, delta, t_value[i]);
        }
        
        //concat(vector_x_query[i], grid_curve_query[i]);
        if (method == 0){
            vector_x_query[i].coords.resize(grid_curve_query[i].coords.size()/2);
            for (int k = 0; k < grid_curve_query[i].coords.size(); k += 2){
                vector_x_query[i].coords[k/2] = grid_curve_query[i].coords[k] + grid_curve_query[i].coords[k+1];
            }
        }
        else{
            vector_x_query[i].coords.resize(grid_curve_query[i].coords.size());
            for (int k = 0; k < grid_curve_query[i].coords.size(); k ++){
                vector_x_query[i].coords[k] = grid_curve_query[i].coords[k];
            }
        }

        vector_x_query[i].hash_id = new long[L];
        hashTables[i]->getBucket(&vector_x_query[i], bucket);
        query_id = vector_x_query[i].hash_id[i];

        for (int j = 0; j < bucket.size(); j++){

            if (method == 0){
                dist = discreteFrechet(*bucket[j]->orig_curve, *query);;
            }
            else{
                dist = continuous_frechet(*bucket[j]->orig_curve, *query);
            }

            //If distance smaller than R, insert it in dR vector
            if (dist <= R){

                //If element already in vector, skip
                if (find(dR.begin(), dR.end(), bucket[j]) != dR.end()){
                    continue;
                }

                dR.push_back(bucket[j]->orig_curve);
            }

        }

    }

    for (int i = 0; i < numOfHashTables; i++){
        delete[] vector_x_query[i].hash_id;
    }

    rR.RDist = dR;
}

//Not as good results as sum of x+y
void LSH_Frechet::concat(pointType& concat, pointType curve){

    for (int k = 0; k < curve.coords.size(); k += 2){

        //concat x, y -> new value k to be stored in vector_x[i][j]
        // -> Doesn't work as good as the sum of x + y

        int int_1, int_2, dec_1, dec_2;
        int_1 = int(curve.coords[k]);
        int_2 = int(curve.coords[k+1]);

        double temp_dec_1 = curve.coords[k] - int_1;
        while (temp_dec_1 / 1.0 != temp_dec_1){
            temp_dec_1 *= 10;
        }
        dec_1 = int(temp_dec_1);

        double temp_dec_2 = curve.coords[k+1] - int_2;
        while (temp_dec_2 / 1.0 != temp_dec_2){
            temp_dec_2 *= 10;
        }
        dec_2 = int(temp_dec_2);

        string s_1 = to_string(int_1), s_2 = to_string(int_2), s_3 = to_string(dec_1), s_4 = to_string(dec_2);
        string int_part = s_1 + s_2;
        string dec_part = s_3 + s_4;

        long int_p = stol(int_part);
        double dec_p = stod(dec_part);

        while (dec_p > 0){
            dec_p /= 10;
        }

        concat.coords[k/2] = int_p + dec_p;
    }

}