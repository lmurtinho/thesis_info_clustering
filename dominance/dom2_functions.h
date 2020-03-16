#ifndef DOM2_FUNCTIONS
#define DOM2_FUNCTIONS

#include "../helper/helper_functions.h"
#include "dominance_heap.h"

int *get_dominant_components(double *data, int n, int dim) {
  int i, j;
  double val, max_val;
  int *ans = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    max_val = data[i*dim];
    ans[i] = 0;
    for (j = 1; j < dim; j++) {
      val = data[i*dim + j];
      if (val > max_val) {
        ans[i] = j;
        max_val = val;
      }
    }
  }
  return ans;
}

double *get_ratios(double *data, int n, int dim) {
  int i, j;
  double val, max_val, sum_vals;
  double *ans = (double *)malloc(sizeof(double) * n);
  for (i = 0; i < n; i++) {
    val = data[i*dim];
    sum_vals = val;
    max_val = val;
    for (j = 1; j < dim; j++) {
      val = data[i*dim + j];
      sum_vals += val;
      if (val > max_val) {
        max_val = val;
      }
    }
    ans[i] = max_val / sum_vals;
  }
  return ans;
}

int *order_data(int *vinds_per_dom, int *vecs_per_dom, int n, int dim) {
  int i, j;
  int *ans = (int *)malloc(sizeof(int) * n);
  int no_idx = 0;
  for (i = 0; i < dim; i++) {
    for (j = 0; j < vecs_per_dom[i]; j++) {
      ans[no_idx] = vinds_per_dom[i*n + j];
      no_idx++;
    }
  }
  return ans;
}

double *rearrange(double *data, int *new_order, int n, int dim) {
  int i;
  double *ans = (double *)malloc(sizeof(double) * n * dim);
  for (i = 0; i < n; i++) {
    copy_to(&data[new_order[i]*dim], &ans[i*dim], dim);
  }
  return ans;
}

double *get_weighted_entropies(double *data, int n, int dim) {
  int i;
  double *ans = (double *)malloc(sizeof(double) * n);
  for (i = 0; i < n; i++) {
    ans[i] = weighted_entropy(&data[i*dim], dim);
  }
  return ans;
}

#endif
