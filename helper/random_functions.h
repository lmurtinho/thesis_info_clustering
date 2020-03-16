#ifndef RANDOM_FUNCTIONS
#define RANDOM_FUNCTIONS

int random_number(int max_number) {
  return rand() % max_number;
}

int find_next_number(double *probs, int n) {
  double r = (double)rand() / (double)(RAND_MAX);
  double sum_probs = 0;
  int i;

  for (i = 0; i < n; i++) {
    sum_probs += probs[i];
    if (sum_probs > r) {
      break;
    }
  }
  return i;
}

int *rd_clustering(double *data, int n, int k, int dim,
                    double *time, int seed) {
    srand(seed);
    clock_t start = clock();
    int i, *assigned = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++) {
      assigned[i] = random_number(k);
    }
    clock_t end = clock();
    time[0] = (double) (end - start) / CLOCKS_PER_SEC;
    printf("C time for rd: %f seconds\n", time[0]);
    return assigned;
  }

double* create_random_array(int n, int dim, int max) {
    int i;
    double *ans = (double *)malloc(sizeof(double) * n * dim);
    for (i = 0; i < n * dim; i++) {
      ans[i] = rand() % max + 1;
    }
    return ans;
  }

int* create_random_int_array(int n, int dim, int max) {
    int i;
    int *ans = (int *)malloc(sizeof(int) * n * dim);
    for (i = 0; i < n * dim; i++) {
      ans[i] = rand() % max + 1;
    }
    return ans;
  }

int* random_sampling(double *values, int n, int k, bool replace) {
  int i, j, *ans = (int *)malloc(sizeof(int) * k);
  double *probs = (double *)malloc(sizeof(double) * n);
  copy_to(values, probs, n);
  normalize_vector(probs, sum(probs, n), n);
  for (i = 0; i < k; i++) {
    ans[i] = find_next_number(probs, n);
    if (!replace) {
      probs[ans[i]] = 0;
      normalize_vector(probs, sum(probs, n), n);
    }
  }
  free(probs);
  return ans;
}

#endif
