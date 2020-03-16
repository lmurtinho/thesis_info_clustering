#ifndef BREG_INIT_FUNCTIONS
#define BREG_INIT_FUNCTIONS

#include "../helper/helper_functions.h"
#include "../helper/random_functions.h"
#include "../dhillon/dhillon_functions.h"

typedef double (*BregmanDivFunc) (double *v1, double *v2, int dim);

typedef double* (*InitClustersFunc) (double *data,
  BregmanDivFunc div_func, int n, int k, int dim);

double* random_init(double *data, BregmanDivFunc div_func,
    int n, int k, int dim) {

  int i, c;
  double *clusters = (double *)malloc(sizeof(double) * k * dim);

  for (i = 0; i < k; i++) {
    c = random_number(n);
    copy_to(&data[c], &clusters[i], dim);
  }

  return clusters;
}

double* pp_init(double *data, BregmanDivFunc div_func,
    int n, int k, int dim) {

  int i, j, prev;
  double sum_divs, new_div;

  double *clusters = (double *)malloc(sizeof(double) * k * dim);
  double *probs = (double *)malloc(sizeof(double) * n);
  double *divs = (double *)malloc(sizeof(double) * n);
  fill_array(divs, BIG_DOUBLE, n);

  // first cluster defined at random
  int nxt = random_number(n);
  copy_to(&data[nxt*dim], clusters, dim);

  // remaining clusters
  for (i = 1; i < k; i++) {

    // update distances
    prev = i-1;
    for (j = 0; j < n; j++) {
      new_div = div_func(&data[j*dim], &clusters[prev*dim], dim);
      if (new_div < divs[j]) {
        divs[j] = new_div;
      }
    }

    // update probabilities
    copy_to(divs, probs, n);
    sum_divs = sum(divs, n);
    normalize_vector(probs, sum_divs, n);

    // find next centroid
    nxt = find_next_number(probs, n);
    copy_to(&data[nxt*dim], &clusters[i*dim], dim);
  }
  free(probs);
  free(divs);
  return clusters;
}

#endif
