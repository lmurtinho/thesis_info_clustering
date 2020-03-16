#include <math.h>

// Returns the squared Euclidean distance between two vectors of size dim
double div_sqd_dst(double *v1, double *v2, int dim) {

  int i;
  double val, ans = 0;

  for (i = 0; i < dim; i++) {
    val = v1[i] - v2[i];
    val *= val;
    ans += val;
  }

  return ans;
}

// Returns the Kullback-Leibler divergence between v1 and v2 (both of size dim)
double div_kl(double *v1, double *v2, int dim) {

  double *norm_v1 = (double *)malloc(sizeof(double) * dim);
  copy_to(v1, norm_v1, dim);
  double norm_v1_sum = sum(norm_v1, dim);
  normalize_vector(norm_v1, norm_v1_sum, dim);

  double *norm_v2 = (double *)malloc(sizeof(double) * dim);
  copy_to(v2, norm_v2, dim);
  double norm_v2_sum = sum(norm_v2, dim);
  normalize_vector(norm_v2, norm_v2_sum, dim);

  double ans = 0, log;
  int i;
  for (i = 0; i < dim; i++) {
    log = log2(norm_v1[i]) - log2(norm_v2[i]);
    ans += norm_v1[i] * log;
  }

  free(norm_v1);
  free(norm_v2);
  return ans;
}

double div_entropy_loss(double *v1, double *v2, int dim) {

  double *v_sum = (double *)malloc(sizeof(double) * dim);
  copy_to(v1, v_sum, dim);
  add_array(v_sum, v2, dim);

  double ans = weighted_entropy(v_sum, dim) - weighted_entropy(v1, dim) -
    weighted_entropy(v2, dim);
  free(v_sum);
  return ans;
}

double div_gini_loss(double *v1, double *v2, int dim) {
  double *v_sum = (double *)malloc(sizeof(double) * dim);
  copy_to(v1, v_sum, dim);
  add_array(v_sum, v2, dim);
  double ans = weighted_gini(v_sum, dim) - weighted_gini(v1, dim) -
    weighted_gini(v2, dim);
  free(v_sum);
  return ans;
}
