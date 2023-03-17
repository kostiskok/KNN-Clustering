#include <fstream>
#include <iostream>

#include "LSH_Vector.h"
#include "Hypercube.h"
#include "LSH_Frechet.h"
#include "clustering.h"
#include "data_handling.h"
#include "functions.h"

using namespace std;

int main(int argc, char *argv[]){

    string inputFile, configFile, outputFile;
    bool inputHasValue = false, configHasValue = false, outputHasValue = false, assignmentHasValue = false, updateHasValue = false;
    bool complete = false, silhouette = false;
    char update, assignment;

    int K=1, L=3, k_LSH=3, M=10, k_HC = 3, probes = 2;
    double delta=1;

    //Read from command line
    for (int i = 0; i < argc; i++){

        string argv_i = argv[i];
        if (argv_i == "-i"){
            inputHasValue = true;
            inputFile = argv[i+1];
        }
        else if (argv_i == "-c"){
            configHasValue = true;
            configFile = argv[i+1];
        }
        else if (argv_i == "-o"){
            outputHasValue = true;
            outputFile = argv[i+1];
        }
        else if (argv_i == "-assignment"){

            string next_argv = argv[i+1];
            if (next_argv == "Classic"){
                assignment = 'c';
                assignmentHasValue = true;
            }
            else if (next_argv == "LSH"){
                assignment = 'l';
                assignmentHasValue = true;
                update = 'v';
                updateHasValue = true;
            }
            else if (next_argv == "Hypercube"){
                assignment = 'h';
                assignmentHasValue = true;
                update = 'v';
                updateHasValue = true;
            }
            else if (next_argv == "Frechet"){
                assignment = 'f';
                assignmentHasValue = true;
                update = 'f';
                updateHasValue = true;
            }
        }
        else if (argv_i == "-update"){

            string next_argv = argv[i+1];
            string next_argv_ = argv[i+2];
            if (next_argv == "Mean" && next_argv_ == "Frechet"){
                update = 'f';
                updateHasValue = true;
            }
            else if (next_argv == "Mean" && next_argv_ == "Vector"){
                update = 'v';
                updateHasValue = true;
            }
        }
        else if (argv_i == "-complete"){
            complete = true;
        }
        else if (argv_i == "-silhouette"){
            silhouette = true;
        }

    }

    if (!inputHasValue){
        cout << "Enter input file path: ";
        cin >> inputFile;
        inputHasValue = true;
    }

    while (!assignmentHasValue){

        string temp;
        cout << "Enter assignment method: ";
        cin >> temp;

        if (temp == "Classic"){
            assignment = 'c';
            assignmentHasValue = true;
        }
        else if (temp == "LSH"){
            assignment = 'l';
            assignmentHasValue = true;
        }
        else if (temp == "Hypercube"){
            assignment = 'h';
            assignmentHasValue = true;
        }
        else if (temp == "Frechet"){
            assignment = 'f';
            assignmentHasValue = true;
        }

    }

    while (!configHasValue){
        cout << "Enter config file path: ";
        cin >> configFile;
        configHasValue = true;
    }

    while (!updateHasValue){

        string temp;
        cout << "Enter update method: ";
        cin >> temp;

        if (temp == "Mean Frechet"){
            update = 'f';
            updateHasValue = true;
        }
        else if (temp == "Mean Vector"){
            update = 'v';
            updateHasValue = true;
        }
    }

    //get parameters config file
    ifstream data;
    data.open(configFile);    //open file
    while(!data){
        cout<<"Could not open file. Please try again\nEnter your input file name\n-->";
        cin>>configFile;
        data.open(configFile); //open file
    }
    string coords;
    while(data>>coords){
        if(coords=="number_of_clusters:"){
            data>>coords;
            K=stoi(coords);
        }
        else if(coords=="number_of_vector_hash_tables:"){
            data>>coords;
            L=stoi(coords);
        }
        else if(coords=="number_of_vector_hash_functions:"){
            data>>coords;
            k_LSH=stoi(coords);
        }
        else if(coords=="max_number_M_hypercube:"){
            data>>coords;
            M=stoi(coords);
        }
        else if(coords=="number_of_hypercube_dimensions:"){
            data>>coords;
            k_HC=stoi(coords);
        }
        else if(coords=="number_of_probes:"){
            data>>coords;
            probes=stoi(coords);
        }
        else if(coords=="delta:"){
            data>>coords;
            delta=stoi(coords);
        }
        
    }
    data.close();

    //Check if the assignments/updates are compatible
    if (assignment == 'f' && update == 'v'){
        cout << "Assignment LSH Frechet and update Mean Vector not possible" << endl;
        cout << "Will compute Mean Frechet instead" << endl;
        update = 'f';
    }
    if ((assignment == 'l' || assignment == 'h') && update == 'f'){
        if (assignment == 'l'){
            cout << "Assignment LSH Vector ";
        }
        else{
            cout << "Assignment Hypercube ";
        }
        cout << "and update Mean Frechet not possible" << endl;
        cout << "Will compute Mean Vector instead" << endl;
        update = 'v';
    }

    //config file

    vector<pointType> points = get_data(inputFile);

    if (update == 'f')
        add_x_values(points);
    int d = points[0].coords.size();

    clustering* cl;

    if (assignment == 'c'){
        cl = new clustering(points, K, d, update);
    }
    else if (assignment == 'l'){
        cl = new clustering(points, K, L, k_LSH, d);
    }
    else if (assignment == 'h'){
        cl = new clustering(points, K, M, probes, k_HC, d);
    }
    else{
        cl = new clustering(points, K, k_LSH, L, d, delta);
    }

    resultCL rCL;
    cl->cluster(rCL, silhouette);

    ofstream outfile(outputFile);

    outfile << "Algorithm: ";
    outfile << "Assignment_";
    if (assignment == 'c')
        outfile << "Lloyds";
    else if (assignment == 'l')
        outfile << "LSH_Vector";
    else if (assignment == 'h')
        outfile << "Hypercube";
    else
        outfile << "LSH_Frechet";
    outfile << ", Update_";
    if (update == 'v')
        outfile << "Mean_Vector";
    else
        outfile << "Mean_Frechet";
    outfile << endl;

    for (int i = 0; i < K; i++){

        outfile << rCL.centroids[i].id;
        outfile << "{size: " << rCL.points_per_centroid[i].size() << ", centroid: ";

        for (int c = 0; c < rCL.centroids[i].coords.size(); c++){
            outfile << rCL.centroids[i].coords[c] << " ";
        }
        outfile << "}" << endl;

    }
    outfile << "clustering_time: " << rCL.time << endl;

    if (silhouette == true){
        outfile << endl << "Silhouette: [";
        for (int i = 0; i < K; i++){
            outfile << rCL.silhouette[i] << ", ";
        }
        outfile << rCL.silhouette[K] << "]" << endl;
    }

    if (complete == true){
        outfile << endl;
        for (int i = 0; i < K; i++){

            outfile << rCL.centroids[i].id << " {";
            for (int c = 0; c < rCL.centroids[i].coords.size(); c++){
                outfile << rCL.centroids[i].coords[c] << " ";
            }

            for (int j = 0; j < rCL.points_per_centroid[i].size(); j++){
                outfile << ", " << rCL.points_per_centroid[i][j]->id;
            }
            outfile << "}" << endl;

        }

        outfile << endl;
    }
    outfile.close();

    delete cl;
}
