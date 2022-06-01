#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <math.h>
#include <ctime>

double eucledian_distance(const double* point_a, const double* point_b,int dim);
double** fill_arr_with_points(int& size,int dim);
void make_random_ndpoint(double* arr, int dim,int lower_bound, int upper_bound);

int main(int argc, char ** argv){
    srand(time(0));
    int n = 0; // number of points in input
    int dim = 2; // dimensions of data
    int K = 5; // We want to identify K different cluseter in our data

    double** data_points = fill_arr_with_points(n,dim);

    int rank, size; 
    double previousMSE = 0;
    double MSE = std::numeric_limits<double>::max(); //MSE should infinity
    double cluster_points[K][dim];

    int number_of_cluster_points[K];

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );


    if(rank == 0){
        for(int i = 0; i < K; i++){
            int random_idx = rand() % n;
            for(int d = 0; d < dim; d++){
                cluster_points[i][d] = data_points[random_idx][d];
            }
            // make_random_ndpoint(cluster_points[i], 2, 0, 3000);
        } 
    }

    MPI_Bcast(&cluster_points,K*dim,MPI_DOUBLE,0,MPI_COMM_WORLD);

    double cluster_points_t[K][dim];
    int number_of_cluster_points_t[K];
    double MSE_t;

    while(MSE != previousMSE){

        previousMSE = MSE;
        MSE_t = 0;

        for(int j = 0; j < K; j++){
            for(int d = 0; d < dim; d++){
                cluster_points_t[j][d] = 0;
            }
            number_of_cluster_points_t[j] = 0;
        }

        double min_dis;
        int closest_cluster_idx;
        for(int i = rank * (n/size); i < (rank+1)*(n/size); i++){ // seperate into block distribution
            double min_dis = eucledian_distance(cluster_points[0],data_points[i],dim);
            int closest_cluster_idx = 0;
            for(int j = 1; j < K; j++){
                //Calculate which centroid is closest to the point
                double curr_dis = eucledian_distance(cluster_points[j],data_points[i],dim);
                if(curr_dis < min_dis){
                    min_dis = curr_dis;
                    closest_cluster_idx = j;
                }
            }

            for(int dim_i = 0; dim_i < dim; dim_i++){
                 cluster_points_t[closest_cluster_idx][dim_i] += data_points[i][dim_i];
            }
           
            number_of_cluster_points_t[closest_cluster_idx]++;
            MSE_t += min_dis;
        }
        MPI_Barrier(MPI_COMM_WORLD); 
        for(int j = 0; j < K; j++){ // find average point for each cluster
                //Calculate which centroid is closest to the point
                MPI_Allreduce(&number_of_cluster_points_t[j], &number_of_cluster_points[j],1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(&cluster_points_t[j],&cluster_points[j], dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                number_of_cluster_points[j] = std::max(number_of_cluster_points[j],1);// to prevent divide by zero

                for(int d = 0; d < dim; d++){ 
                    cluster_points[j][d] /= number_of_cluster_points[j];
                }
            }

        MPI_Allreduce(&MSE_t, &MSE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    
    if(rank == 0){
        for(int j = 0; j < K; j++){
            for(int d = 0; d < dim; d++){
                std::cout << cluster_points[j][d];
                if(d != dim - 1){
                    std::cout << ",";
                }
            }
            std::cout << "\n";
        }
    }
    
    MPI_Finalize();
    
    
    for(int i = 0; i < n; i++){
         delete [] data_points[i];
    }
    delete [] data_points;
    return 0;
}


double eucledian_distance(const double* point_a, const double* point_b, int dim){
    double distance = 0.0;
    double raw_dis;

    for(int i = 0; i < dim; i++){
        raw_dis = point_a[i] - point_b[i];
        distance += raw_dis * raw_dis;
    }
    
    return sqrt(distance);
}


void make_random_ndpoint(double * arr,int dim, int lower_bound, int upper_bound){
    for(int i = 0; i < dim; i++){
        arr[i] = rand() % upper_bound;
    }
}


double** fill_arr_with_points(int& size,int dim){

    std::ifstream myFile("./test/blobs.txt");
    static double** arr = new double*[1000000];
    for(int i = 0; i < 1000000; i++){
        arr[i] = new double[dim];
    }
    std::string line;
    double val;
    if(myFile.good())
    {
        size = 0;
        while(std::getline(myFile, line)){
            std::stringstream stream(line);
            for(int i = 0; i < dim - 1; i++){
                stream >> arr[size][i];
                stream.ignore(); 
            }
            stream >> arr[size][dim - 1];
            size++;       
        }

    }
    
    myFile.close();
    return arr;
}
