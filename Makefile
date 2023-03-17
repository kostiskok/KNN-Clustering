CC = g++
CFLAGS = -g -Wall -I. -I./fred_files

OBJ_CL = main_Cluster.o clustering.o LSH_Frechet.o LSH_Vector.o Hypercube.o HashTable_LSH.o HashTable_Hypercube.o functions.o data_handling.o continuous_frechet.o
OBJ_LSH_F = main_Search.o LSH_Frechet.o LSH_Vector.o Hypercube.o HashTable_LSH.o HashTable_Hypercube.o functions.o data_handling.o continuous_frechet.o

OBJ_FRED = fred_files/config.o fred_files/curve.o fred_files/frechet.o fred_files/interval.o fred_files/point.o fred_files/simplification.o

cluster: $(OBJ_CL) $(OBJ_FRED)
	$(CC) -o $@ $^ $(CFLAGS)

search: $(OBJ_LSH_F) $(OBJ_FRED)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o fred_files/*.o lsh cube cluster search