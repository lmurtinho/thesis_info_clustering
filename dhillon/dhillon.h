#ifndef DHILLON
#define DHILLON

#include "dhillon_functions.h"

typedef bool (*InitClusterFunc) (double *data, double *clusters,
  double *norm_data, double *data_logs, int n, int k, int dim);

bool expectation(double *norm_data, double *data_logs,
                  int *assigned, double *clusters,
                  int n, int k, int dim) {
  int i, j, curr, new;
  double dist, d, c_sum;

  double* cluster_logs = get_logs(clusters, k, dim);
  int *vecs_per_cluster = (int *)malloc(sizeof(int)*k);
  fill_array_int(vecs_per_cluster, 0, k);

  bool changed = false;                      // no change made so far

  for (i = 0; i < n; i++) {              // for each vector
    dist = BIG_DOUBLE;                 // current distance is largest possible
    curr = assigned[i];                 // get vector's current assignment
    new = -1;                            // new assignment is non-existent
    // dist = kl_div(norm_data, data_logs, cluster_logs, i, new, dim);
    for (j = 0; j < k; j++) {            // for each cluster
      // CHANGE LATER
      d = kl_div(norm_data, data_logs, cluster_logs, i, j, dim);  // get distance between cluster and vector
      if (d < dist)  {                    // if smaller than smallest distance
        new = j;                         // new assignment is current cluster
        vecs_per_cluster[j]++;
        dist = d;                        // new distance is current distance
      }
      if (d >= BIG_DOUBLE) {
        printf("infinite distance!\n");
      }
    }
    if (new != curr) {                   // if assignment changed
      changed = true;                   // set change to true
      assigned[i] = new;                // change assignment to new cluster
    }
    if (new == -1) {
      printf("not assigned!\n");
    }
  }

  for (i = 0; i < k; i++) {
    if (vecs_per_cluster[i] == 0) {
      j = (int) (rand() % n);
      assigned[j] = i;
      changed = true;
    }
  }

  free(vecs_per_cluster);
  free(cluster_logs);
  return changed;
}

// K-MEANS WITH DHILLON INITIALIZATION
int* di_clustering(double *data, int n, int k, int dim,
                    int max_iter, double *t, int seed) {

  srand(seed);
  clock_t start = clock();
  int iter = 0;
  int *assigned = (int *)malloc(sizeof(int) * n);
  int i, j, c;
  double *clusters = (double *)malloc(sizeof(double) * k * dim), pe;

  // assign initial centers (Dhillon implementation)
  initial_di(data, clusters, assigned, n, k, dim);

  double *norm_data = (double *)malloc(sizeof(double) * n * dim);
  copy_to(data, norm_data, n * dim);
  normalize_array(norm_data, n, dim);

  double *data_logs = get_logs_from_normal(norm_data, n * dim);

  bool changed = true;

  while(changed && (iter < max_iter) ) {
    iter++;
    changed = expectation(norm_data, data_logs, assigned, clusters, n, k, dim);
    maximization(data, assigned, clusters, n, k, dim);
  }
  free(clusters);
  free(norm_data);
  free(data_logs);

  clock_t end = clock();
  t[0] = (double) (end - start) / CLOCKS_PER_SEC;
  printf("C time for di: %f seconds\n", t[0]);
  return assigned;
}

#endif
