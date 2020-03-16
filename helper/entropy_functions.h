#ifndef ENTROPY_FUNCTIONS
#define ENTROPY_FUNCTIONS

#include <math.h>

// Returns the weighted Gini of dim-sized double array v
double weighted_gini(double *v, int dim) {
  double ans = 0, v_sum = sum(v, dim);
  if (v_sum != 0) {
    int i;
    for (i = 0; i < dim; i++) {
      ans += v[i] * (1 - (v[i] / v_sum));
    }
  }
  return ans;
}

double weighted_entropy(double *v, int dim) {
  double ans = 0, v_sum = sum(v, dim);
  if (v_sum != 0) {
    int i;
    for (i = 0; i < dim; i++) {
      if (v[i] != 0) {
        ans -= v[i] * log2(v[i] / v_sum);
      }
    }
  }
  return ans;
}

double* clusters_entropy(double *clusters, int k, int dim) {
  int i;

  double* results = (double *)malloc(sizeof(double) * k);
  for (i = 0; i < k; i++) {
    results[i] = weighted_entropy(&clusters[i*dim], dim);
  }

  return results;
}

double partition_entropy(double *clusters, int k, int dim)
  {
    double* results = clusters_entropy(clusters, k, dim);
    double ans = sum(results, k);
    free(results);
    return ans;
  }

#endif
