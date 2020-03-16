#ifndef STRONG_CORESETS
#define STRONG_CORESETS
#include <math.h>

double sqd_eucl_dist(double *v1, double *v2, int dim) {
  double ans = 0, temp;
  int i;
  for (i = 0; i < dim; i++) {
    temp = v1[i] - v2[i];
    temp *= temp;
    ans += temp;
  }
  return ans;
}

// initial centers found according to ++ implementation
void initial_pp(double *data, double *clusters, double *distances,
    int *assigned, int n, int k, int dim, double factor) {
  double d, sum_distances;
  int curr_cluster, i, j;

  double *probs = (double *)malloc(sizeof(double) * n);

  fill_array(distances, (double)INFINITY, n);
  fill_array_int(assigned, 0, n);

  curr_cluster = rand() % n;
  copy_to(&data[curr_cluster * dim], clusters, dim);

  for (i = 1; i < k; i++) {
    for (j = 0; j < n; j++) {
      d = sqd_eucl_dist(&data[j * dim], &clusters[(i-1)*dim], dim);
      d *= factor;
      if (d < distances[j]) {
        distances[j] = d;
        assigned[j] = i - 1;
      }
    }

    copy_to(distances, probs, n);
    sum_distances = sum(distances, n);
    normalize_vector(probs, sum_distances, n);

    curr_cluster = find_next_number(probs, n);

    copy_to(&data[curr_cluster * dim], &clusters[i * dim], dim);
  }

  i = k-1;
  for (j = 0; j < n; j++) {
    d = sqd_eucl_dist(&data[j * dim], &clusters[i * dim], dim);
    d *= factor;
    if (d < distances[j]) {
      distances[j] = d;
      assigned[j] = i;
    }
  }

  free(probs);
  return;
}

int* get_vecs_per_cluster(int *assigned, int n, int k) {
  int i, *ans = (int *)malloc(sizeof(int) * k);
  fill_array_int(ans, 0, n);
  for (i = 0; i < n; i++) {
    ans[assigned[i]]++;
  }
  return ans;
}

double *get_norm_data(double *data, int n, int dim) {
  double *ans = (double *)malloc(sizeof(double) * n * dim);
  copy_to(data, ans, n * dim);
  normalize_array(ans, n, dim);
  return ans;
}

double* get_distances_per_cluster(double *distances, int *assigned,
  int n, int k) {
  int i;
  double *ans = (double *)malloc(sizeof(double) * k);
  fill_array(ans, 0, k);
  for (i = 0; i < n; i++) {
    ans[assigned[i]] += distances[i];
  }
  return ans;
}

double* get_sensitivities(double *data, double *norm_data, double *data_logs,
    double *clusters, int n, int k, int dim, double factor) {

  int i, c;
  double v, alpha = 16 * (log2(k) + 2), // Step 1
    *ans = (double *)malloc(sizeof(double) * n),
    *cluster_logs = (double *)malloc(sizeof(double) * k * dim),
    *distances = (double *)malloc(sizeof(double) * n);

  int *assigned = (int *)malloc(sizeof(int) * n);
  initial_pp(data, clusters, distances, assigned, n, k, dim, factor);

  double *distances_per_cluster = get_distances_per_cluster(distances, assigned,
      n, k);

  double distance_avg = sum(distances, n) / n;

  for (i = 0; i < n; i ++) {
    c = assigned[i];
    v = alpha * distances[i] / distance_avg;
    v += 2 * alpha * distances_per_cluster[c] / (k * distance_avg);
    v += (4 * n) / k;
    ans[i] = v;
  }

  free(assigned);
  free(cluster_logs);
  free(distances);
  free(distances_per_cluster);

  return ans;
}

int* coreset_construction(double *data, double *norm_data, double *data_logs,
    double *clusters, int n, int k, int dim, int m, double factor) {

  double *sensitivities = get_sensitivities(data, norm_data, data_logs,
    clusters, n, k, dim, factor);

  int *ans = random_sampling(sensitivities, n, m, false);

  return ans;
}

void kmeans(double *data, double *norm_data, double *data_logs,
    double *clusters, int n, int k, int dim, int max_iter, int *iter,
    int *iters, double *cluster_retrieval, double *times_to_cluster,
    clock_t initial_time, double factor) {

  int i, j, c, *assigned = (int *)malloc(sizeof(int) * n), nxt = 0;
  clock_t end;
  fill_array_int(assigned, -1, n);
  iter[0] = 0;
  double *cluster_logs = (double *)malloc(sizeof(double) * k * dim),
    *distances = (double *)malloc(sizeof(double) * n);

  initial_pp(data, clusters, distances, assigned, n, k, dim, factor);

  // store initial clusters
  if (iters[nxt] == 0) {
    end = clock();
    copy_to(clusters, &cluster_retrieval[nxt * k * dim], k * dim);
    times_to_cluster[nxt] = (double) (end - initial_time) / CLOCKS_PER_SEC;
    nxt++;
  }

  bool changed = true;

  while(changed && (iter[0] < max_iter) ) {
    iter[0]++;
    changed = expectation(norm_data, data_logs, assigned, clusters, n, k, dim);
    maximization(data, assigned, clusters, n, k, dim);

    // store clusters in array
    if (iters[nxt] == iter[0]) {
      end = clock();
      copy_to(clusters, &cluster_retrieval[nxt * k * dim], k * dim);
      times_to_cluster[nxt] = (double) (end - initial_time) / CLOCKS_PER_SEC;
      nxt++;
    }
  }

  // STORE FINAL CLUSTERS IN ARRAY
  end = clock();
  iters[nxt] = iter[0];
  copy_to(clusters, &cluster_retrieval[nxt * k * dim], k * dim);
  times_to_cluster[nxt] = (double) (end - initial_time) / CLOCKS_PER_SEC;

  free(cluster_logs);
  free(distances);
}

int* co_clustering(double *data, int n, int k, int dim, int m,
  int max_iter, int *iter, char *data_name, int exp) {
  clock_t initial_time = clock(), start = initial_time, end;
  int iters[] = {0, 1, 5, 10, 100}, n_iters = 5;
  double *norm_data = get_norm_data(data, n, dim),
         *data_logs = get_logs_from_normal(norm_data, n * dim),
         *clusters  = (double *)malloc(sizeof(double) * k * dim),
         *cluster_retrieval = (double *)malloc(sizeof(double) * k * dim * n_iters),
         *times_to_cluster = (double *)malloc(sizeof(double) * n_iters);
  int *assigned = (int *)malloc(sizeof(int) * n);
  double full_time;
  int i;
  double lambda = find_min(data, n * dim);
  double factor = 1 / (2 * lambda);

int *coreset_indices = coreset_construction(data, norm_data, data_logs,
    clusters, n, k, dim, m, factor);

  double *coreset_data = get_indices(data, coreset_indices, n, dim, m);
  double *coreset_norm = get_indices(norm_data, coreset_indices, n, dim, m);
  double *coreset_logs = get_indices(data_logs, coreset_indices, n, dim, m);

  kmeans(coreset_data, coreset_norm, coreset_logs, clusters, m, k, dim,
    max_iter, iter, iters, cluster_retrieval, times_to_cluster, initial_time,
    factor);

  FILE *results = fopen("results.csv", "a");

  for (i = 0; i < n_iters; i++) {
    start = clock();
    fill_array_int(assigned, -1, n);
    copy_to(&cluster_retrieval[i * k * dim], clusters, k * dim);
    expectation(norm_data, data_logs, assigned, clusters, n, k, dim);
    end = clock();
    full_time = times_to_cluster[i] + (double) (end - start) / CLOCKS_PER_SEC;
    maximization(data, assigned, clusters, n, k, dim);
    fprintf(results, "%d,%s,coresets,%d,%d,%d,%f,%f\n",
      exp, data_name, k, m, iters[i],
      partition_entropy(clusters, k, dim), full_time);
    if (iters[i] == iter[0]) {
      break;
    }
  }
  fclose(results);

  free(norm_data);
  free(data_logs);
  free(clusters);
  free(coreset_indices);
  free(coreset_data);
  free(coreset_norm);
  free(coreset_logs);
  free(cluster_retrieval);

  clock_t end_time = clock();
  printf("Full time: %f seconds\n", (double) (end_time - initial_time) / CLOCKS_PER_SEC);
  return assigned;
}

#endif
