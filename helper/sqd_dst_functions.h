#ifndef SQD_DST
#define SQD_DST

// Returns the sum of the squared elements in a double array
double squared_sum(double *v, int dim) {
  double ans = 0;
  int i;
  for (i = 0; i < dim; i++) {
    ans += (v[i] * v[i]);
  }
  return ans;
}

double* get_centroids(double *clusters, int *vecs_per_cluster,
  int k, int dim) {
  double *centroids = (double *)malloc(sizeof(double) * k * dim);
  copy_to(clusters, centroids, k * dim);
  int i, j, n;
  for (i = 0; i < k; i++) {
    n = vecs_per_cluster[i];
    if (n > 0) {
      for (j = 0; j < dim; j++) {
        centroids[i*dim + j] /= n;
      }
    }
  }
  return centroids;
}

double sqd_dst(double *u, double *v, int dim) {
  double ans = 0, val;
  int i;
  for (i = 0; i < dim; i++) {
    val = u[i] - v[i];
    val *= val;
    ans += val;
  }
  return ans;
}

double cluster_sqd_dst(double *data, int *assigned, int n, int c, int dim) {
  int i;
  double *sum_vec = (double *)malloc(sizeof(double)*dim);
  fill_array(sum_vec, 0, dim);
  double sqd_sum = 0;
  double sum_vec_sqd = 0;
  int size = 0;
  double ans;

  for (i = 0; i < n; i++) {
    if (assigned[i] == c) {
      sqd_sum += squared_sum(&data[i*dim], dim);
      add_array(sum_vec, &data[i*dim], dim);
      size++;
    }
  }

  if (size > 0) {
    sum_vec_sqd = squared_sum(sum_vec, dim);
    ans = sqd_sum - (sum_vec_sqd / size);
  }

  free(sum_vec);
  return ans;
}

double* clusters_sqd_dst(double *data, int *assigned, int n, int k, int dim) {

  int i;

  double *results = (double *)malloc(sizeof(double) * k);

  for (i = 0; i < k; i++) {
    results[i] = cluster_sqd_dst(data, assigned, n, i, dim);
  }

  return results;

}

double partition_sqd_dst(double *data, int *assigned, int n, int k, int dim)
  {
    double* results = clusters_sqd_dst(data, assigned, n, k, dim);
    double ans = sum(results, k);
    free(results);
    return ans;
  }

#endif
