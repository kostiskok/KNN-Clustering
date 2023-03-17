#include <iostream>
#include <string>
#include <chrono>
#include <map>
#include <algorithm>
#include <random>

#include "clustering.h"
#include "functions.h"

using namespace std;

//clustering: lloyd's algorithm
clustering::clustering(vector<pointType>& p, int K, int d, char update_method): 
    points(p), K(K), d(d), update(update_method), centroids(K), silhouette(K+1, 0){

    method = 'c';

    initialize_common();
    
}

//clustering: range search with LSH
clustering::clustering(vector<pointType>& p, int K, int L, int k, int d): 
    points(p), K(K), L(L), k(k), d(d), centroids(K), silhouette(K+1, 0){

    method = 'l';
    update = 'v';

    lsh_v = new LSH_Vector(points, k, L, d, 6, L2);

    initialize_common();

    for (int i = 0; i < points.size(); i++){
        isAssigned.insert(pair<string, bool>(points[i].id, false));
    }

}

//clustering: range search with Hypercube
clustering::clustering(vector<pointType>& p, int K, int M, int probes, int k, int d): 
    points(p), K(K), M(M), probes(probes), k(k), d(d), centroids(K), silhouette(K+1, 0){

    method = 'h';
    update = 'v';

    hc = new Hypercube(points, k, d, 6, M, probes, L2);

    initialize_common();

    for (int i = 0; i < points.size(); i++){
        isAssigned.insert(pair<string, bool>(points[i].id, false));
    }

}

//clustering: range search with Frechet
clustering::clustering(vector<pointType>& c, int K, int k, int L, int d_p, double delta): 
    points(c), K(K), L(L), k(k), d(d_p), centroids(K), silhouette(K+1, 0){

    method = 'f';
    update = 'f';

    lsh_f = new LSH_Frechet(points, k, L, 6, delta, "discrete");
    d = points[0].coords.size();

    initialize_common();

    for (int i = 0; i < points.size(); i++){
        isAssigned.insert(pair<string, bool>(points[i].id, false));
    }

}

//Initialize all data structrures common for all three algorithms
void clustering::initialize_common(){

    for (int i = 0; i < K; i++){
        centroids[i].id = ("CLUSTER-" + to_string(i+1));
        centroids[i].coords.resize(d);
        for (int j = 0; j < d; j++){
            centroids[i].coords[j] = 0;
        }
    }

    points_per_centroid.resize(K);

    temp_cluster.coords.resize(d);

}

clustering::~clustering(){

    if (method == 'l'){
        delete lsh_v;

        /*for (int i = 0; i < K; i++){
            delete[] centroids[i].hash_id;
        }*/

    }
    else if (method == 'h'){
        delete hc;
    }
    else if (method == 'f'){
        delete lsh_f;

        //????
    }
    

}

void clustering::cluster(resultCL& rCL, bool sil){

    bool flag = false;

    auto start = chrono::high_resolution_clock::now();

    k_means(); //for each algorithm, find initial centroids with k_means++

    while(1){

        //Initialize per loop
        for (int i = 0; i < K; i++){
            points_per_centroid[i].clear();
        }

        //assignment
        if (method == 'c'){
            assignment_lloyds();
        }
        else if (method == 'l'){
            assignment_rangeSearch_LSH_Vector();
        }
        else if (method == 'h'){
            assignment_rangeSearch_Hypercube();
        }
        else{
            assignment_rangeSearch_LSH_Frechet();
        }

        limit = getMinDistance(centroids, update);

        //update
        if (update == 'v'){
            flag = update_Vector();
            if (flag == false){
                break;
            }
        }
        else{
            flag = update_Frechet();
            if (flag == false){
                break;
            }
        }

    }

    auto stop=chrono::high_resolution_clock::now();
    auto duration=chrono::duration<double>(stop-start);

    if (sil){

        //calculate silhouettes
        double average_stotal = 0;
        for (int i = 0; i < K; i++){

            double average_s = 0;
            double a_i, b_i, temp;
            for (int j = 0; j < points_per_centroid[i].size(); j++){

                if (update == 'v')
                    a_i = L2(&centroids[i], points_per_centroid[i][j], d);
                else
                    a_i = discreteFrechet(centroids[i], *points_per_centroid[i][j]);

                b_i = numeric_limits<double>::max();

                for (int other_c = 0; other_c < K; other_c++){

                    if (other_c == i) continue;

                        if (update == 'v'){
                            temp = L2(&centroids[other_c], points_per_centroid[i][j], d);
                        }
                        else{
                            temp = discreteFrechet(centroids[other_c], *points_per_centroid[i][j]);
                        }
                    

                    if (temp < b_i)
                        b_i = temp;

                }
                
                if (a_i < b_i)
                    temp = 1 - a_i/b_i;
                else if(a_i == b_i)
                    temp = 0;
                else
                    temp = b_i/a_i - 1;

                average_s += temp;

            }

            average_s /= points_per_centroid[i].size();
            silhouette[i] = average_s; //average per cluster
            average_stotal += average_s;

        }
        average_stotal /= K;
        silhouette[K] = average_stotal; //average of all clusters

    }

    rCL.centroids = centroids;
    rCL.points_per_centroid = points_per_centroid;
    rCL.time = duration.count();

    if (sil){
        rCL.silhouette = silhouette;
    }

}

void clustering::assignment_lloyds(){

    int index;
    double distance;
    
    pointType result;
    result.coords.resize(d);

    //Assign each point to their closer centroid

    for(int i=0;i<points.size();i++){

        double min = numeric_limits<double>::max(), dist;
        int index;
        for (int j = 0; j < K; j++){

            if (update == 'v')
                dist = L2(&points[i], &centroids[j], d);
            else
                dist = discreteFrechet(points[i], centroids[j]);


            if (dist < min){
                min = dist;
                index = j;
            }

        }

        points_per_centroid[index].push_back(&points[i]);
    }

}

void clustering::assignment_rangeSearch_LSH_Vector(){
    vector<vector<pointType*>> temp;
    resultR rR;
    int count;

    pointType temp_cluster;
    temp_cluster.coords.resize(d);

    //isAssigned(point_name) = true if point.id == point_name is assigned to some cluster

    //initializations for each loop
    for (int i = 0; i < points.size(); i++){
        isAssigned[points[i].id] = false;
    }

    //assignment -> reverse approach:
    // - find the points closer to each centroid and assign them to it
    // - LSH::rangeSearch(radius) (if a point lies in > 1 balls, then compare to each centroid and assign it)
    // - multiply radius by 2, till no new elements are assigned to new balls (approximate)
    // - for each point that hasn't been assigned, compare with each centroid

    //start with radius = min(dist between centers)/2
    double radius = numeric_limits<double>::max();
    
    for (int i = 0; i < K; i++){
        for (int j = i+1; j < K; j++){

            double distance = L2(&centroids[i], &centroids[j], d);
            if (radius > distance){
                radius = distance;
            }

        }
    }
    radius /= 2;

    while(1){

        count = 0;
        temp.clear(); //in each iteration clear the temp vector

        //For each centroid, find the radius-close points to it and store them temporary in the temp vector
        // (temp[i] stores all the points close to centroid[i])
        for (int i = 0; i < K; i++){

            lsh_v->rangeSearch(&centroids[i], radius, rR);

            temp.push_back(rR.RDist);

        }

        for (int i = 0; i < K; i++){ //for each vector stored in temp
            for (int j = 0; j < temp[i].size(); j++){ //for each point stored in temp[i]

                if (isAssigned[temp[i][j]->id] == false){

                    //Search if the unassigned point is concurrently at two or more Rs, and
                    // assign it to the closest
                    int centroid_index = i;
                    double distance = L2(temp[i][j], &centroids[i], d);
                    
                    for (int other_c = 0; other_c < K; other_c++){

                        if (other_c == i) continue;
                        if (find(temp[other_c].begin(), temp[other_c].end(), temp[i][j]) != temp[other_c].end()){

                            double temp_distance = L2(temp[i][j], &centroids[other_c], d);
                            if (temp_distance < distance){
                                distance = temp_distance;
                                centroid_index = other_c;
                            }

                        }

                    }

                    points_per_centroid[centroid_index].push_back(temp[i][j]);
                    isAssigned[temp[i][j]->id] = true;
                    count++;

                }

            }
        }

        //if at max one new point is inserted in the centroid, stop
        if (count <= 5)
            break;

        radius *= 2;

    }

    //For each unassigned point, find closest centroid
    for (int i = 0; i < points.size(); i++){

        if (isAssigned[points[i].id] == false){

            double min = numeric_limits<double>::max();
            int centroid_index;

            for (int j = 0; j < K; j++){

                double distance = L2(&points[i], &centroids[j], d);

                if (min > distance){
                    min = distance;
                    centroid_index = j;
                }

            }

            points_per_centroid[centroid_index].push_back(&points[i]);
            isAssigned[points[i].id] == true;

        }

    }

}

void clustering::assignment_rangeSearch_Hypercube(){

    vector<vector<pointType*>> temp;
    resultR rR;
    int count;

    pointType temp_cluster;
    temp_cluster.coords.resize(d);

    //isAssigned(point_name) = true if point.id == point_name is assigned to some cluster

    //initializations for each loop
    for (int i = 0; i < points.size(); i++){
        isAssigned[points[i].id] = false;
    }

    //assignment -> reverse approach:
    // - find the points closer to each centroid and assign them to it
    // - LSH::rangeSearch(radius) (if a point lies in > 1 balls, then compare to each centroid and assign it)
    // - multiply radius by 2, till no new elements are assigned to new balls (approximate)
    // - for each point that hasn't been assigned, compare with each centroid

    //start with radius = min(dist between centers)/2
    double radius = numeric_limits<double>::max();
    for (int i = 0; i < K; i++){
        for (int j = i+1; j < K; j++){

            double distance = L2(&centroids[i], &centroids[j], d);
            if (radius > distance){
                radius = distance;
            }

        }
    }
    radius /= 2;

    while(1){

        count = 0;
        temp.clear(); //in each iteration clear the temp vector

        //For each centroid, find the radius-close points to it and store them temporary in the temp vector
        // (temp[i] stores all the points close to centroid[i])
        for (int i = 0; i < K; i++){

            hc->rangeSearch(&centroids[i], radius, rR);

            temp.push_back(rR.RDist);

        }

        for (int i = 0; i < K; i++){ //for each vector stored in temp
            for (int j = 0; j < temp[i].size(); j++){ //for each point stored in temp[i]

                if (isAssigned[temp[i][j]->id] == false){

                    //Search if the unassigned point is concurrently at two or more Rs, and
                    // assign it to the closest
                    int centroid_index = i;
                    double distance = L2(temp[i][j], &centroids[i], d);
                    
                    for (int other_c = 0; other_c < K; other_c++){

                        if (other_c == i) continue;
                        if (find(temp[other_c].begin(), temp[other_c].end(), temp[i][j]) != temp[other_c].end()){

                            double temp_distance = L2(temp[i][j], &centroids[other_c], d);
                            if (temp_distance < distance){
                                distance = temp_distance;
                                centroid_index = other_c;
                            }

                        }

                    }

                    points_per_centroid[centroid_index].push_back(temp[i][j]);
                    isAssigned[temp[i][j]->id] = true;
                    count++;

                }

            }
        }

        //if at max one new point is inserted in the centroid, stop
        if (count <= 5)
            break;

        radius *= 2;

    }

    //For each unassigned point, find closest centroid
    for (int i = 0; i < points.size(); i++){

        if (isAssigned[points[i].id] == false){

            double min = numeric_limits<double>::max();
            int centroid_index;

            for (int j = 0; j < K; j++){

                double distance = L2(&points[i], &centroids[j], d);

                if (min > distance){
                    min = distance;
                    centroid_index = j;
                }

            }

            points_per_centroid[centroid_index].push_back(&points[i]);
            isAssigned[points[i].id] == true;

        }

    }


}

void clustering::assignment_rangeSearch_LSH_Frechet(){
    vector<vector<pointType*>> temp;
    resultR rR;
    int count;

    pointType temp_cluster;
    temp_cluster.coords.resize(d);

    //isAssigned(point_name) = true if point.id == point_name is assigned to some cluster

    //initializations for each loop
    for (int i = 0; i < points.size(); i++){
        isAssigned[points[i].id] = false;
    }

    //assignment -> reverse approach:
    // - find the points closer to each centroid and assign them to it
    // - LSH::rangeSearch(radius) (if a point lies in > 1 balls, then compare to each centroid and assign it)
    // - multiply radius by 2, till no new elements are assigned to new balls (approximate)
    // - for each point that hasn't been assigned, compare with each centroid

    //start with radius = min(dist between centers)/2
    double radius = numeric_limits<double>::max();

    for (int i = 0; i < K; i++){
        for (int j = i+1; j < K; j++){

            double distance = discreteFrechet(centroids[i], centroids[j]);
            if (radius > distance){
                radius = distance;
            }

        }
    }
    radius /= 2;

    while(1){

        count = 0;
        temp.clear(); //in each iteration clear the temp vector

        //For each centroid, find the radius-close points to it and store them temporary in the temp vector
        // (temp[i] stores all the points close to centroid[i])
        for (int i = 0; i < K; i++){

            lsh_f->rangeSearch(&centroids[i], radius, rR);

            temp.push_back(rR.RDist);

        }

        for (int i = 0; i < K; i++){ //for each vector stored in temp
            for (int j = 0; j < temp[i].size(); j++){ //for each point stored in temp[i]

                if (isAssigned[temp[i][j]->id] == false){

                    //Search if the unassigned point is concurrently at two or more Rs, and
                    // assign it to the closest
                    int centroid_index = i;
                    double distance = discreteFrechet(*temp[i][j], centroids[i]);
                    
                    for (int other_c = 0; other_c < K; other_c++){

                        if (other_c == i) continue;
                        if (find(temp[other_c].begin(), temp[other_c].end(), temp[i][j]) != temp[other_c].end()){

                            double temp_distance = discreteFrechet(*temp[i][j], centroids[other_c]);
                            if (temp_distance < distance){
                                distance = temp_distance;
                                centroid_index = other_c;
                            }

                        }

                    }

                    points_per_centroid[centroid_index].push_back(temp[i][j]);
                    isAssigned[temp[i][j]->id] = true;
                    count++;

                }

            }
        }

        //if at max one new point is inserted in the centroid, stop
        if (count <= 5)
            break;

        radius *= 2;

    }

    //For each unassigned point, find closest centroid
    for (int i = 0; i < points.size(); i++){

        if (isAssigned[points[i].id] == false){

            double min = numeric_limits<double>::max();
            int centroid_index;

            for (int j = 0; j < K; j++){

                double distance = discreteFrechet(points[i], centroids[j]);

                if (min > distance){
                    min = distance;
                    centroid_index = j;
                }

            }

            points_per_centroid[centroid_index].push_back(&points[i]);
            isAssigned[points[i].id] == true;

        }

    }

}

bool clustering::update_Vector(){
    // - for each cluster, calculate mean point value and assign it as the new cluster
    // - if there is little change in the centroids, stop

    bool flag = false;

    for (int i = 0; i < K; i++){

        for (int c = 0; c < d; c++){
            temp_cluster.coords[c] = 0;
        }

        //For each point.id stored in each cluster:
        for (int j = 0; j < points_per_centroid[i].size(); j++){

            //sum all coords..
            for (int c = 0; c < d; c++){
                temp_cluster.coords[c] += points_per_centroid[i][j]->coords[c];
            }

        }
        //..so as to calculate the mid value (aka the new centroid)
        for (int k = 0; k < d; k++){
            temp_cluster.coords[k] /= points_per_centroid[i].size();
        }

        //if the new centroid is far from the old one ( > 0.1), there will be a new loop
        if (L2(&temp_cluster, &centroids[i], d) > limit){
            flag = true;
        }
        
        assignPoint(centroids[i], temp_cluster);

    }

    return flag;
}

bool clustering::update_Frechet(){

    random_device r;
    default_random_engine gen(r());
    bool flag = false;
    vector<pointType*> tempPoints_per_centroid;
    vector<pointType*> meanTemp;

    for (int i = 0; i < K; i++){
        
        tempPoints_per_centroid.clear();

        for (int j = 0; j < points_per_centroid[i].size(); j++){

            tempPoints_per_centroid.push_back(points_per_centroid[i][j]);

        }

        while (tempPoints_per_centroid.size() > 1){

            while (tempPoints_per_centroid.size() > 1){
                //get two random curves

                uniform_int_distribution<int> uni_i(0, tempPoints_per_centroid.size()-1);

                int first, second;
                first = uni_i(gen);

                do {
                    second = uni_i(gen);
                } while (second == first);

                //find the mean curve of these 2
                for (int k = 0; k < d; k++){
                    tempPoints_per_centroid[first]->coords[k] = (tempPoints_per_centroid[first]->coords[k] + tempPoints_per_centroid[second]->coords[k])/2;
                }
                meanTemp.push_back(tempPoints_per_centroid[first]);

                if (first > second){
                    tempPoints_per_centroid.erase(tempPoints_per_centroid.begin() + first);
                    tempPoints_per_centroid.erase(tempPoints_per_centroid.begin() + second);
                }
                else{
                    tempPoints_per_centroid.erase(tempPoints_per_centroid.begin() + second);
                    tempPoints_per_centroid.erase(tempPoints_per_centroid.begin() + first);
                }

            }

            //tempPoints_per_centroid.resize(meanTemp.size());
            for (int j = 0; j < meanTemp.size(); j++){
                tempPoints_per_centroid.push_back(meanTemp[j]);
            }

            meanTemp.clear();

        }

        if (discreteFrechet(*tempPoints_per_centroid[0], centroids[i]) > limit){
            flag = true;
        }

        assignPoint(centroids[i], *tempPoints_per_centroid[0]);

    }

    return flag;

}

//k_means++ to find initial centroids
void clustering::k_means(){

    int nearestCentroid, index;
    double max, prob_sum;
    long double sum;

    bool* isCentroid;
    isCentroid=new bool[points.size()];
    for(int i=0;i<points.size();i++){
        isCentroid[i]=false;
    }

    vector<double> dist(points.size());
    vector<double> probability(points.size());
    
    random_device r;
    default_random_engine gen(r());
    uniform_int_distribution<int> uni_i(0, points.size()-1);
    int first_centroid = uni_i(gen);

    //First centroid is completely at random
    assignPoint(centroids[0],points[first_centroid]);
    isCentroid[first_centroid]=true;

    //For the rest of centroids:
    for(int j=1;j<K;j++){
        prob_sum=0;
        sum=0;
        for(int i=0;i<points.size();i++){ //calculate each point's minimum distance from every assigned centroid
            if(isCentroid[i]==false){
                dist[i]=minimum_distance(points[i],nearestCentroid,j);
            }
        }
        //and then calculate each point's probability to become a new centroid
        max=*max_element(std::begin(dist), std::end(dist));
        for(int i=0;i<points.size();i++){
            if(isCentroid[i]==false){
                sum+=pow(dist[i],2)/max;
            }
        }
        for(int i=0;i<points.size();i++){ //NOTE: Probability is too small + it doesn't add up to 1 for all points
            if(isCentroid[i]==false){
                probability[i]=dist[i]/sum;
                prob_sum+=probability[i];
            }
        }

        //new weighted random centroid
        uniform_real_distribution<double> uni_d(0,prob_sum);
        double random=uni_d(gen);
        for(int i=0;i<points.size();i++){
            if(isCentroid[i]==false){ 
                if(random<=probability[i]){
                    assignPoint(centroids[j],points[i]);
                    isCentroid[i]=true;
                    break;
                }
                random-=probability[i]; 
            }
        }

    }
    delete[] isCentroid;
}

//point target has now the coords of the point source
void clustering::assignPoint(pointType& target, pointType source){
    for(int i=0;i<target.coords.size();i++){
        target.coords[i]=source.coords[i];
    }
}

//returns p's minimum distance from the first |limit| centroids
double clustering::minimum_distance(pointType p, int& index, int limit){
    double min=numeric_limits<double>::max();
    double distance;
    for(int i=0;i<limit;i++){

        if (update == 'v')
            distance=L2(&p, &centroids[i],d);
        else
            distance=discreteFrechet(p, centroids[i]);

        if(distance<min){
            min=distance;
            index=i;
        }
    }
    return min;
}