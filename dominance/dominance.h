#ifndef DOM
#define DOM

#include "dominance_functions.h"
#include "dominance_heap.h"

double *transform_data(double *data, int n, int k, int dim) {
  double *ans = (double *)malloc(sizeof(double)*n*k),
    *data_sum = array_sum(data, n, dim), rest, val;
  int* le = largest_elements(data_sum, dim, k-1), i, j, p;

  for (i = 0; i < n; i++) {
    rest = sum(&data[i*dim], dim);

    for (j = 0; j < k-1; j++) {
      p = le[j];
      val = data[i*dim + p];
      rest -= val;
      ans[i*k+j] = val;
    }
    ans[i*k+(k-1)] = rest;
  }

  return ans;
}

int* rg_main(double *data, int n, int k, int dim) {
    int i, j, d;

    int *assigned = (int *)malloc(sizeof(int) * n);

    // FIND CLUSTER ASSIGNMENTS FOR K = DIM AND VECTOR RATIOS
    double *clusters = (double *)malloc(sizeof(double) * k * dim);
    double *ratio = (double *)malloc(sizeof(double) * n);
    dom_first_pass(data, clusters, assigned, ratio, n, k, dim);

    // SORT INDICES OF CLUSTERS BY ENTROPY
    double *original_results = clusters_entropy(clusters, dim, dim);
    double o_sum = sum(original_results, dim);
    normalize_vector(original_results, o_sum, dim);
    int *c_indices = largest_elements(original_results, dim, dim);

    // SORT INDICES OF VECTORS BY ENTROPY
    int *v_indices = largest_elements(ratio, n, n);

    // PARTITION VECTOR INDICES BY FIRST ASSIGNMENT CLUSTER
    int *vinds_per_cluster = (int *)malloc(sizeof(int) * dim * n);
    int *vecs_per_cluster = (int *)malloc(sizeof(int) * dim);
    first_partition(vinds_per_cluster, vecs_per_cluster, v_indices, assigned,
      n, dim);

    // FIND ADDITIONAL CLUSTERS PER COMPONENT
    int *additional_clusters = assign_clusters(c_indices,
      original_results, vecs_per_cluster, k, dim);

    double *cluster_data = (double *)malloc(sizeof(double)*n*dim);
    int *cluster_assignment = (int *)malloc(sizeof(int)*n);
    int first_cluster = 0;

    for (i = 0; i < dim; i++) {
      d = c_indices[i];
      if (vecs_per_cluster[d] > 0) {
        fill_array_int(cluster_assignment, -1, n);
        get_vecs_for_cluster(data, cluster_data, &vinds_per_cluster[d * n], dim,
          vecs_per_cluster[d]);
        cluster_assign_ext(cluster_data, cluster_assignment, first_cluster,
          additional_clusters[d], vecs_per_cluster[d], dim);

        for (j = 0; j < vecs_per_cluster[d]; j++) {
          assigned[vinds_per_cluster[d*n + j]] = cluster_assignment[j];
        }
        first_cluster += additional_clusters[d];
      }
    }

    free(clusters);
    free(ratio);
    free(original_results);
    free(c_indices);
    free(v_indices);
    free(vinds_per_cluster);
    free(vecs_per_cluster);
    free(additional_clusters);
    free(cluster_data);
    free(cluster_assignment);
    return assigned;
}

int* rg_clustering(double *data, int n, int k,
                    int dim, double *time) {

  clock_t start = clock();

  int *assigned = (int *)malloc(sizeof(int)*n);

  if (k < dim) {
    double *new_data = transform_data(data, n, k, dim);
    assigned = rg_main(new_data, n, k, k);
    free(new_data);
  }
  else {
    assigned = rg_main(data, n, k, dim);
  }

  clock_t end = clock();

  time[0] = (double) (end - start) / CLOCKS_PER_SEC;
  printf("C time for rg: %f seconds\n", time[0]);
  return assigned;
}

#endif
