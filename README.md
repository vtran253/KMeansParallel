# CPSC479-Kmeans-Parallel
**Project 2:**
A parallel implementation of the k-means problem using Lloyd's heuristics. Implemented with Message passing interface (MPI) 

## Contributers
Cody Mangham cmangham1@csu.fullerton.edu

Vivian Tran vtran2535@csu.fullerton.edu

Carson Carpenter carson.carpenter7@csu.fullerton.edu

## Runnning
Ensure you have MPI installed on machine
1. `make`
2. `mpirun -np <process num> ./main `
This will print the coordinates of the centroids in the terminal
## Testing
A python script is included to visualize the data set and thhe cluster points
### Checking
1. Run Program
2. Copy points into the cluster_points.txt 
3. Run `python3 ./test/create_dataset.py`
