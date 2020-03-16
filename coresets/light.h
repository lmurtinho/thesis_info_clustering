#ifndef LIGHT_CORESET
#define LIGHT_CORESET

#include "../helper/helper_functions.h"
#include "../dhillon/dhillon.h"
#include "coresets.h"

void get_logs2 (double *a, double *logs, int size) {
  int i;
  for (i = 0; i < size; i++) {
    logs[i] = log2(a[i]);
  }
}

double opt_kl_div(double *norm1, double *log1, double *log2, int dim) {
  double ans = 0, log;
  int i;

  for (i = 0; i < dim; i++) {
    log = log1[i] - log2[i];
    ans += norm1[i] * log;
  }
  return ans;
}

double simple_kl_div(double *v1, double *v2, int dim) {

  double *n1 = (double *)malloc(sizeof(double) * dim);
  double *n2 = (double *)malloc(sizeof(double) * dim);

  copy_to(v1, n1, dim);
  copy_to(v2, n2, dim);

  double s1 = sum(n1, dim);
  normalize_vector(n1, s1, dim);
  double s2 = sum(n2, dim);
  normalize_vector(n2, s2, dim);


  double *l1 = get_logs_from_normal(n1, dim);
  double *l2 = get_logs_from_normal(n2, dim);

  double ans = 0, log;
  int i;

  for (i = 0; i < dim; i++) {
    log = l1[i] - l2[i];
    ans += n1[i] * log;
  }

  free(n1);
  free(n2);
  free(l1);
  free(l2);
  return ans;
}

// initial centers found according to ++ implementation
void initial_pp_kl(double *data, double *clusters, double *norm_data,
    double *data_logs, int n, int k, int dim) {

  double *probs         = (double *)malloc(sizeof(double) * n);
  double *norm_clusters = (double *)malloc(sizeof(double) * k * n);
  double *distances = (double *)malloc(sizeof(double) * n);
  double *cluster_logs = (double *)malloc(sizeof(double) * k * dim);

  double d, sum_distances;
  int curr_cluster, i, j;

  fill_array(distances, (double)INFINITY, n);

  curr_cluster = rand() % n;
  copy_to(&data[curr_cluster * dim], clusters, dim);
  copy_to(clusters, norm_clusters, dim);
  normalize_array(norm_clusters, 1, dim);

  for (i = 0; i < dim; i++) {
    cluster_logs[i] = log2(norm_clusters[i]);
  }

  for (i = 1; i < k; i++) {
    for (j = 0; j < n; j++) {
      d = kl_div(norm_data, data_logs, cluster_logs, j, i-1, dim);
      d *= d;
      if (d < distances[j]) {
        distances[j] = d;
      }
    }

    copy_to(distances, probs, n);
    sum_distances = sum(distances, n);
    normalize_vector(probs, sum_distances, n);

    curr_cluster = find_next_number(probs, n);

    copy_to(&data[curr_cluster * dim], &clusters[i * dim], dim);
    copy_to(&clusters[i * dim], &norm_clusters[i * dim], dim);
    normalize_array(&norm_clusters[i * dim], 1, dim);

    for (j = i*dim; j < (i*dim + dim); j++) {
      cluster_logs[j] = log2(norm_clusters[j]);
    }

  }
  free(probs);
  free(norm_clusters);
  free(distances);
  free(cluster_logs);
  return;
}

CORESET_SAMPLE light_coreset (double *data, int n, int k, int m, int dim,
    int seed) {
  int i;

  srand(seed);
  double *mean = (double *)malloc(sizeof(double) * dim);
  double *sqd_dst = (double *)malloc(sizeof(double) * n);
  double *qs = (double *)malloc(sizeof(double) * n);
  double *all_weights = (double *)malloc(sizeof(double) * n);
  double *weights = (double *)malloc(sizeof(double) * m);
  double *coreset = (double *)malloc(sizeof(double) * m * dim);

  double *norm_data = (double *)malloc(sizeof(double) * n * dim);
  copy_to(data, norm_data, n * dim);
  normalize_array(norm_data, n, dim);

  double *logs_norm = (double *)malloc(sizeof(double) * n * dim);
  copy_to(norm_data, logs_norm, n * dim);
  logalize_vector(logs_norm, n * dim);

  fill_array(mean, 0, dim);
  for (i = 0; i < n; i++) {
    add_array(mean, &data[i*dim], dim);
  }
  for (i = 0; i < dim; i++) {
    mean[i] /= n;
  }

  double *logs_mean = get_logs(mean, 1, dim);

  double c = 1 / (2 * n);
  for (i = 0; i < n; i++) {
    sqd_dst[i] = opt_kl_div(&norm_data[i*dim], &logs_norm[i*dim],
                            logs_mean, dim);
    sqd_dst[i] *= sqd_dst[i];
    if (i == 0) {
    }
  }
  double sum_sqd_dst = sum(sqd_dst, n);
  sum_sqd_dst /= 2;

  for (i = 0; i < n; i++) {
    qs[i] = c + sqd_dst[i] / sum_sqd_dst;
  }

  for (i = 0; i < n; i++) {
    all_weights[i] = 1 / (m * qs[i]);
  }

  int *samples = random_sampling(qs, n, m, false);

  for (i = 0; i < m; i++) {
    weights[i] = all_weights[samples[i]];
  }
  double weight_sum = sum(weights, m);
  normalize_vector(weights, weight_sum, m);

  for(i = 0; i < m; i++) {
    copy_to(&data[samples[i]*dim], &coreset[i*dim], dim);
  }

  CORESET_SAMPLE ans;
  ans.samples = coreset;
  ans.weights = weights;

  free(norm_data);
  free(logs_norm);
  free(mean);
  free(logs_mean);
  free(sqd_dst);
  free(qs);
  free(all_weights);
  free(samples);

  return ans;
}

double* bregman_kl_pp(double *data, double *weights,
                int n, int k, int dim, int max_iter) {

  // create arrays for cluster assignment and cluster vectors
  int *assigned = (int *)malloc(sizeof(int) * n);
  double *clusters = (double *)malloc(sizeof(double) * k * dim), pe;
  int i, j, c, iter;

  double *norm_data = (double *)malloc(sizeof(double) * n * dim);
  copy_to(data, norm_data, n * dim);
  normalize_array(norm_data, n, dim);

  double *data_logs = get_logs_from_normal(norm_data, n * dim);

  initial_pp_kl(data, clusters, norm_data, data_logs, n, k, dim);

  bool changed = expectation(norm_data, data_logs, assigned,
                              clusters, n, k, dim);
  iter = 0;

  while(changed && (iter < max_iter) ) {
    iter++;
    weighted_maximization(data, assigned, clusters, weights, n, k, dim);
    changed = expectation(norm_data, data_logs, assigned,
                            clusters, n, k, dim);
  }
  free(assigned);
  free(norm_data);
  free(data_logs);
  return clusters;
}


int *lc_clustering (double *data, int n, int k,
                    int dim, int m, int max_iter, double *t, int seed) {

  srand(seed);
  clock_t start = clock();

  int *assigned = (int *)malloc(sizeof(int) * n);

  CORESET_SAMPLE coreset = light_coreset(data, n, k, m, dim, seed);

  double *clusters = bregman_kl_pp(coreset.samples, coreset.weights,
                                      m, k, dim, max_iter);

  double *norm_data = (double *)malloc(sizeof(double) * n * dim);
  copy_to(data, norm_data, n * dim);
  normalize_array(norm_data, n, dim);

  double *data_logs = get_logs_from_normal(norm_data, n * dim);

  bool changed = expectation(norm_data, data_logs, assigned, clusters,
                             n, k, dim);
  free(clusters);
  free(norm_data);
  free(data_logs);
  free(coreset.samples);
  free(coreset.weights);

  clock_t end = clock();
  t[0] = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("C time for lc: %f seconds\n", t[0]);
  return assigned;
}

#endif
