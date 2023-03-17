# KNN-Clustering
Find K-nearest neighbors or k Clusters from a given timeseries

## Input format
```
$./search –i <input file> –q <query file> –k <int> -L <int> -M <int> -probes <int> -ο <output file> -algorithm <LSH or Hypercube or Frechet> -metric <discrete or continuous | only for –algorithm Frechet> -delta <double>

$./cluster –i <input file> –c <configuration file> -o <output file> -update <Mean Frechet or Mean Vector> –assignment <Classic or LSH or Hypercube or LSH_Frechet> -complete <optional> -silhouette <optional>
```
