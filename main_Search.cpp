#include <fstream>
#include <iostream>

#include "LSH_Vector.h"
#include "Hypercube.h"
#include "LSH_Frechet.h"
#include "data_handling.h"
#include "functions.h"

using namespace std;

int main(int argc, char *argv[]){

    int k = 4, L = 5, delta = 1, M = 10, probes = 2;
    char algorithm, metric;
    string inputFile, queryFile, outputFile;
    bool inputHasValue = false, queryHasValue = false, outputHasValue = false, algorithmHasValue = false, metricHasValue = false;;

    LSH_Vector* LSH_V;
    Hypercube* HC;
    LSH_Frechet* LSH_F;

    //Read from command line
    for (int i = 0; i < argc; i++){

        string argv_i = argv[i];
        if (argv_i == "-i"){
            inputHasValue = true;
            inputFile = argv[i+1];
        }
        else if (argv_i == "-q"){
            queryHasValue = true;
            queryFile = argv[i+1];
        }
        else if (argv_i == "-k"){
            k = atoi(argv[i+1]);
        }
        else if (argv_i == "-L"){
            L = atoi(argv[i+1]);
        }
        else if (argv_i == "-probes"){
            probes = atoi(argv[i+1]);
        }
        else if (argv_i == "-o"){
            outputHasValue = true;
            outputFile = argv[i+1];
        }
        else if (argv_i == "-algorithm"){

            string next_argv = argv[i+1];
            if (next_argv == "LSH"){
                algorithm = 'l';
                algorithmHasValue = true;
            }
            else if (next_argv == "Hypercube"){
                algorithm = 'h';
                algorithmHasValue = true;
            }
            else if (next_argv == "Frechet"){
                algorithm = 'f';
                algorithmHasValue = true;
            }
        }
        else if (argv_i == "-metric"){

            string next_argv = argv[i+1];
            if (next_argv == "discrete"){
                metric = 'd';
                metricHasValue = true;
            }
            else if (next_argv == "continuous"){
                metric = 'c';
                metricHasValue = true;
            }
        }
        else if (argv_i == "-delta"){
            string delta_string = argv[i+1];
            delta = stod(delta_string);
        }

    }

    if (!inputHasValue){
        cout << "Enter input file path: ";
        cin >> inputFile;
        inputHasValue = true;
    }

    while (!algorithmHasValue){
        string temp;
        cout << "Enter algorithm: ";
        cin >> temp;
        if (temp == "LSH"){
            algorithm = 'l';
            algorithmHasValue = true;
        }
        else if (temp == "Hypercube"){
            algorithm = 'h';
            algorithmHasValue = true;
        }
        else if (temp == "Frechet"){
            algorithm = 'f';
            algorithmHasValue = true;
        }
    }

    while (algorithm == 'f' && !metricHasValue){
        string temp;
        cout << "Enter metric for -algorithm Frechet: ";
        cin >> temp;
        if (temp == "discrete"){
            metric = 'd';
            metricHasValue = true;
        }
        else if (temp == "continuous"){
            metric = 'c';
            metricHasValue = true;
        }

    }

    vector<pointType> curves = get_data(inputFile);
    if ((algorithm == 'f') && (metric == 'd')){
        add_x_values(curves);
    }

    int d = curves[0].coords.size();

    if (algorithm == 'l'){
        LSH_V = new LSH_Vector(curves, k, L, d, 6, L2);
    }
    else if (algorithm == 'h'){
        HC = new Hypercube(curves, k, d, 6, M, probes, L2);
    }
    else{
        if (metric == 'd'){
            LSH_F = new LSH_Frechet(curves, k, L, 6, delta, "discrete");
        }
        else{
            LSH_F = new LSH_Frechet(curves, k, L, 6, delta, "continuous");
        }
    }

    

    while(1){

        if (!queryHasValue){
            cout << "Enter path for query file: ";
            cin >> queryFile;
            queryHasValue = true;
        }

        if (!outputHasValue){
            cout << "Enter path for output file: ";
            cin >> outputFile;
            outputHasValue = true;
        }

        ofstream outfile(outputFile);
        vector<pointType> queries = get_data(queryFile);

        if ((algorithm == 'f') && (metric == 'd')){
            add_x_values(queries);
        }

        double tApproximateAverage = 0;
        double tTrueAverage = 0;
        double maf = -1;

        for (int q = 0; q < queries.size(); q++){

            cout << "Searching for nearest curve for query " << queries[q].id << "..." << endl;

            resultNN rNN;
            if (algorithm == 'l'){
                LSH_V->nearestNeighbours(&queries[q], 1, rNN);
            }
            else if (algorithm == 'h'){
                HC->nearestNeighbours(&queries[q], 1, rNN);
            }
            else{
                LSH_F->nearestNeighbours(&queries[q], 1, rNN);
            }

            outfile << "Query: " << queries[q].id << endl;
            if (algorithm == 'l'){
                outfile << "Algorithm: LSH_Vector";
            }
            else if(algorithm == 'h'){
                outfile << "Algorithm: Hypercube";
            }
            else{
                if (metric == 'd')
                    outfile << "Algorithm: LSH_Frechet_Discrete";
                else
                    outfile << "Algorithm: LSH_Frecher_Continuous";
            }
            outfile << endl;
            outfile << "Approximate Nearest neighbor: " << rNN.Dist[0].second << endl;
            outfile << "True Nearest Neighbor: " << rNN.trueDist[0].second << endl;
            outfile << "distanceApproximate: " << rNN.Dist[0].first << endl;
            outfile << "distanceTrue: " << rNN.trueDist[0].first << endl;

            tApproximateAverage += rNN.tDist;
            tTrueAverage += rNN.tTrue;

            double current_af = rNN.Dist[0].first / rNN.trueDist[0].first;
            if (current_af > maf){
                maf = current_af;
            }

            outfile << endl;

            cout << "Done" << endl;

        }

        tApproximateAverage /= queries.size();
        tTrueAverage /= queries.size();

        outfile << "tApproximateAverage: " << tApproximateAverage << endl;
        outfile << "tTrueAverage: " << tTrueAverage << endl;
        outfile << "MAF: " << maf << endl;

        outfile.close();

        break;

    }

    if (algorithm == 'l'){
        delete LSH_V;
    }
    else if (algorithm == 'h'){
        delete HC;
    }
    else{
        //!!! CONTINUOUS
        delete LSH_F;
    }

    
}