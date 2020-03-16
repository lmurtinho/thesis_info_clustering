#ifndef DHILLON_FUNCTIONS
#define DHILLON_FUNCTIONS

#include "../helper/helper_functions.h"

// find cluster of vector v in dominance algorithm
int get_cluster(double *v, int *le, int dim, int k) {
  double v_sum = sum(v, dim);

  double max = -INFINITY;
  int idx, i, max_idx = -1;

  for (i = 0; i < k-1; i++) {
    idx = le[i];
    v_sum -= v[idx];
    if (v[idx] > max) {
      max = v[idx];
      max_idx = i;
    }
  }

  if (v_sum > max) {
    max_idx = k-1;
  }

  return max_idx;
}

void dominance_small_k(double *data, int *assigned, int n, int k, int dim) {
  int i;
  int l = (k > dim)? dim: k;

  double* data_sum = array_sum(data, n, dim);
  int* le = largest_elements(data_sum, dim, l);

  for (i = 0; i < n; i++) {
    assigned[i] = get_cluster(&data[i * dim], le, dim, l);
  }

  free(data_sum);
  free(le);
}

double kl_div(double *norm_data, double *data_logs, double *cluster_logs,
              int n, int k, int dim) {
  double ans = 0, log;
  int i;
  int v = n * dim;
  int c = k * dim;

  for (i = 0; i < dim; i++) {
    log = data_logs[v + i] - cluster_logs[c + i];
    ans += norm_data[v + i] * log;
  }

  return ans;
}

void initial_di(double *data, double *clusters, int *assigned,
                     int n, int k, int dim) {
  int iter[] = {0};
  dominance_small_k(data, assigned, n, k, dim);

  if (k > dim) {

    int *clusters_per_component = (int *)malloc(sizeof(int)*dim);
    fill_array_int(clusters_per_component, k / dim, dim);
    int addtl_clusters = k - int_sum(clusters_per_component, dim);
    int i = 0;

    while (addtl_clusters > 0) {
      clusters_per_component[i]++;
      i++;
      addtl_clusters--;
    }

    int *initial_cluster = (int *)malloc(sizeof(int) * dim);
    int *current_cluster = (int *)malloc(sizeof(int) * dim);
    initial_cluster[0] = 0;
    current_cluster[0] = 0;

    for (i = 1; i < dim; i++) {
      initial_cluster[i] = initial_cluster[i-1] + clusters_per_component[i-1];
      current_cluster[i] = initial_cluster[i];
    }

    int c;

    for (i = 0; i < n; i++) {
      c = assigned[i];
      assigned[i] = current_cluster[c];
      current_cluster[c]++;
      if (current_cluster[c] >=
          (initial_cluster[c] + clusters_per_component[c]) ) {
        current_cluster[c] = initial_cluster[c];
      }
    }

    free(clusters_per_component);
    free(initial_cluster);
    free(current_cluster);
  }

  maximization(data, assigned, clusters, n, k, dim);
}

#endif
