#ifndef LIGHT_CORESET
#define LIGHT_CORESET

#include "../helper/helper_functions.h"
#include "../dhillon/dhillon.h"
#include "coresets.h"
#include "../helper/sqd_dst_functions.h"

// Weighted maximization: recalculate centers according to vectors assignment
void km_weighted_maximization(double *data, int *assigned, double *centroids,
                           double *weights, int n, int k, int dim) {
  int i, j, c;
  double w;
  fill_array(centroids, 0, k*dim);
  double *c_weights = (double *)malloc(sizeof(double) * k);
  fill_array(c_weights, 0, k);
  // for each vector
  for (i = 0; i < n; i++) {
    c = assigned[i];
    w = weights[i];
    c_weights[c] += w;
    for (j = 0; j < dim; j++) {
      centroids[c * dim + j] += (data[i * dim + j] * w);
    }
  }

  // divide by weights
  for (i = 0; i < k; i++) {
    for (j = 0; j < dim; j++) {
      centroids[i*dim + j] /= c_weights[i];
    }
  }
  return;
}

bool bregman_expectation(double *data, double *centroids, int *assigned,
    int n, int k, int dim) {

  bool changed = false;
  int i, j, curr, new;
  double min_div, div;

  for (i = 0; i < n; i++) {
    min_div = BIG_DOUBLE;
    curr = assigned[i];
    new = -1;
    for (j = 0; j < k; j++) {
      div = sqd_eucl_dist(&data[i*dim], &centroids[j*dim], dim);
      if (div < min_div) {
        new = j;
        min_div = div;
      }
      else if (div >= BIG_DOUBLE) {
        printf("infinite distance\n");
      }
    }
    if (new != curr) {
      changed = true;
      assigned[i] = new;
    }
    if (new == -1) {
      printf("not assigned\n");
    }
  }
  return changed;
}

// initial centers found according to ++ implementation
void lkm_initial_pp_kl(double *data, double *centroids,
    int n, int k, int dim) {

  double *probs     = (double *)malloc(sizeof(double) * n);
  double *distances = (double *)malloc(sizeof(double) * n);

  double d, sum_distances;
  int c, i, j;

  fill_array(distances, (double)INFINITY, n);

  c = rand() % n;
  copy_to(&data[c * dim], centroids, dim);

  for (i = 1; i < k; i++) {
    for (j = 0; j < n; j++) {
      d = sqd_eucl_dist(&data[j*dim], &centroids[(i-1)*dim], dim);
      if (d < distances[j]) {
        distances[j] = d;
      }
    }

    copy_to(distances, probs, n);
    sum_distances = sum(distances, n);
    normalize_vector(probs, sum_distances, n);

    c = find_next_number(probs, n);

    copy_to(&data[c * dim], &centroids[i * dim], dim);
    }
  free(probs);
  free(distances);
}

CORESET_SAMPLE lkm_light_coreset (double *data, int n, int m, int dim,
    int seed) {
  int i;

  srand(seed);

  double *distances = (double *)malloc(sizeof(double) * n);
  double *qs = (double *)malloc(sizeof(double) * n);
  double *all_weights = (double *)malloc(sizeof(double) * n);
  double *weights = (double *)malloc(sizeof(double) * m);
  double *coreset = (double *)malloc(sizeof(double) * m * dim);

  double *mean = mean_vector(data, n, dim);

  double c = 1 / (2 * n);
  for (i = 0; i < n; i++) {
    distances[i] = sqd_eucl_dist(&data[i*dim], mean, dim);
  }
  double sum_dist = sum(distances, n);
  sum_dist /= 2;

  for (i = 0; i < n; i++) {
    qs[i] = c + (distances[i] / sum_dist);
    all_weights[i] = 1 / (m * qs[i]);
  }

  int *samples = random_sampling(qs, n, m, false);

  for (i = 0; i < m; i++) {
    weights[i] = all_weights[samples[i]];
  }

  for(i = 0; i < m; i++) {
    copy_to(&data[samples[i]*dim], &coreset[i*dim], dim);
  }

  CORESET_SAMPLE ans;
  ans.samples = coreset;
  ans.weights = weights;

  free(mean);
  free(distances);
  free(qs);
  free(all_weights);
  free(samples);

  return ans;
}

double* lkm_bregman_kl_pp(double *data, double *weights,
                int n, int k, int dim, int max_iter) {

  // create arrays for cluster assignment and cluster vectors
  int *assigned = (int *)malloc(sizeof(int) * n);
  double *clusters = (double *)malloc(sizeof(double) * k * dim), pe;
  int i, j, c, iter;

  lkm_initial_pp_kl(data, clusters, n, k, dim);

  bool changed = bregman_expectation(data, clusters, assigned,
                                      n, k, dim);
  iter = 0;
  while(changed && (iter < max_iter) ) {
    iter++;
    km_weighted_maximization(data, assigned, clusters, weights, n, k, dim);
    changed = bregman_expectation(data, clusters, assigned,
                                  n, k, dim);
  }
  free(assigned);
  return clusters;
}


int *lkm_lc_clustering (double *data, int n, int k,
                    int dim, int m, int max_iter, double *t, int seed) {

  srand(seed);
  clock_t start = clock();

  int *assigned = (int *)malloc(sizeof(int) * n);

  CORESET_SAMPLE coreset = lkm_light_coreset(data, n, m, dim, seed);
  double *centroids = lkm_bregman_kl_pp(coreset.samples, coreset.weights,
                                      m, k, dim, max_iter);
  bool changed = bregman_expectation(data, centroids, assigned,
                                      n, k, dim);
  free(centroids);
  free(coreset.samples);
  free(coreset.weights);

  clock_t end = clock();
  t[0] = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("C time for lc: %f seconds\n", t[0]);
  return assigned;
}


#endif
