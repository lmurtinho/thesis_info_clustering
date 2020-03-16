#ifndef BREGMAN
#define BREGMAN

#include "../dhillon/dhillon_functions.h"
#include "../dhillon/dhillon.h"
#include "../coresets/light.h"

double* pp_init(double *data, int n, int k, int dim) {

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
      new_div = simple_kl_div(&data[j*dim], &clusters[prev*dim], dim);
      if (new_div < divs[j]) {
        divs[j] = new_div;
      }
    }

    // update probabilities
    copy_to(divs, probs, n);
    sum_divs = sum(divs, n);
    // printf("loss for %d clusters: %.20e\n", prev, sum_divs);
    normalize_vector(probs, sum_divs, n);

    // find next centroid
    nxt = find_next_number(probs, n);
    copy_to(&data[nxt*dim], &clusters[i*dim], dim);
  }
  free(probs);
  free(divs);
  return clusters;
}

// K-MEANS WITH ++ INITIALIZATION
int* bh_clustering(double *data, int n, int k, int dim,
                    int max_iter, double *t, int seed) {

  srand(seed);
  clock_t start = clock();
  int iter = 0;
  int *assigned = (int *)malloc(sizeof(int) * n);
  int i, j, c;
  double *clusters = pp_init(data, n, k, dim);

  double *norm_data = (double *)malloc(sizeof(double) * n * dim);
  copy_to(data, norm_data, n * dim);
  normalize_array(norm_data, n, dim);

  double *data_logs = get_logs_from_normal(norm_data, n * dim);

  bool changed = expectation(norm_data, data_logs, assigned, clusters, n, k, dim);

  while(changed && (iter < max_iter) ) {
    iter++;
    maximization(data, assigned, clusters, n, k, dim);
    changed = expectation(norm_data, data_logs, assigned, clusters, n, k, dim);
  }
  free(clusters);
  free(norm_data);
  free(data_logs);

  clock_t end = clock();
  t[0] = (double) (end - start) / CLOCKS_PER_SEC;
  printf("C time for bh: %f seconds\n", t[0]);
  return assigned;
}

#endif
