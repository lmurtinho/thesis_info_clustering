#ifndef HELPER_FUNCTIONS
#define HELPER_FUNCTIONS

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#define BIG_DOUBLE (INFINITY)

typedef struct CORESET_SAMPLE {
  double *samples;
  double *weights;
} CORESET_SAMPLE;

typedef struct HEAP {
  int id;
  int ini;
  int end;
  double val;
  double *ini_vec;
  double *end_vec;
  double *heap_vec;
  int prev;
  int next;
} HEAP;

int *index_array(int n) {
  int i;
  int *ans = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    ans[i] = i;
  }
  return ans;
}

double find_min (double *a, int n) {
  int i;
  double ans = BIG_DOUBLE;
  for (i = 0; i < n; i++) {
    if (a[i] < ans) {
      ans = a[i];
    }
  }
  return ans;
}

double *cum_sum(double *a, int dim) {
  int i;
  double *ans = (double *)malloc(sizeof(double) * dim);
  ans[0] = a[0];
  for (i = 1; i < dim; i++) {
    ans[i] = ans[i-1] + a[i];
  }
  return ans;
}

int *cum_sum_int(int *a, int dim) {
  int i;
  int *ans = (int *)malloc(sizeof(int) * dim);
  ans[0] = a[0];
  for (i = 1; i < dim; i++) {
    ans[i] = ans[i-1] + a[i];
  }
  return ans;
}

// Returns the sum of elements in an array
double sum(double *v, int dim) {
  double ans = 0;
  int i;
  for (i = 0; i < dim; i++) {
    ans += v[i];
  }
  return ans;
}

int int_sum(int *v, int dim) {
  int ans = 0, i;
  for (i = 0; i < dim; i++) {
    ans += v[i];
  }
  return ans;
}

// Calculates sum of all columns in a "2d" array
double* array_sum(double *a, int n, int dim) {
  int i, j;
  double* ans = (double *)malloc(sizeof(double) * dim);
  for (i = 0; i < dim; i++) {
    ans[i] = 0;
    for (j = 0; j < n; j++) {
      ans[i] += a[j*dim + i];
    }
  }
  return ans;
}

void fill_array(double *a, double v, int dim) {
  int i;
  for (i = 0; i < dim; i++) {
    a[i] = v;
  }
}

void fill_array_int(int *a, int v, int dim) {
  int i;
  for (i = 0; i < dim; i++) {
    a[i] = v;
  }
}

void add_array(double *dest, double *add, int dim) {
  int i;
  for (i = 0; i < dim; i++) {
    dest[i] += add[i];
  }
}

double *mean_vector(double *a, int n, int dim) {
  double *mean = (double *)malloc(sizeof(double) * dim);
  fill_array(mean, 0, dim);
  int i;
  for (i = 0;  i < n; i++) {
    add_array(mean, &a[i*dim], dim);
  }
  for (i = 0; i < dim; i++) {
    mean[i] /= n;
  }
  return mean;
}

// Maximization phase: recalculate centers according to vectors assignment
void maximization(double *data, int *assigned, double *clusters,
                 int n, int k, int dim) {
  int i, j, l, cluster;
  fill_array(clusters, 0, k*dim);
  for (i = 0; i < n; i++) {
    cluster = assigned[i];
    for (j = 0; j < dim; j++) {
      clusters[cluster * dim + j] += data[i * dim + j];
    }
  }
  return;
}

// Weighted maximization: recalculate centers according to vectors assignment
void weighted_maximization(double *data, int *assigned, double *clusters,
                           double *weights, int n, int k, int dim) {
  int i, j, l, cluster;
  fill_array(clusters, 0, k*dim);
  // for each vector
  for (i = 0; i < n; i++) {
    cluster = assigned[i];
    for (j = 0; j < dim; j++) {
      clusters[cluster * dim + j] += weights[i] * data[i * dim + j];
    }
  }
  return;
}

void logalize_vector(double *v, int n) {
  int i;
  for (i = 0; i < n; i++) {
    v[i] = log2(v[i]);
  }
}

void normalize_vector(double *v, double v_sum, int dim) {
  int i;
  for (i = 0; i < dim; i++) {
    v[i] /= v_sum;
  }
  return;
}

void normalize_array(double *a, int n, int dim) {
  int i;
  double v_sum;
  for (i = 0; i < n; i++) {
    v_sum = sum(&a[i * dim], dim);
    normalize_vector(&a[i * dim], v_sum, dim);
  }
  return;
}

void copy_to(double *src, double *dst, int dim) {
  int i;
  for (i = 0; i < dim; i++) {
    dst[i] = src[i];
  }
  return;
}

void copy_to_int(int *src, int *dst, int dim) {
  int i;
  for (i = 0; i < dim; i++) {
    dst[i] = src[i];
  }
  return;
}

double* get_logs_from_normal(double *a, int size) {
  int i;
  double *logs_a = (double *)malloc(sizeof(double) * size);
  for (i = 0; i < size; i++) {
    logs_a[i] = log2(a[i]);
  }
  return logs_a;
}


double* get_logs(double *a, int n, int dim) {
  int i;
  double* temp = (double *)malloc(sizeof(double) * n * dim);
  copy_to(a, temp, n * dim);
  normalize_array(temp, n , dim);
  double *logs_a = get_logs_from_normal(temp, n * dim);
  free(temp);
  return logs_a;
}

int compare(double *array, int i, int j) {
  int ans = 0;
  if (array[i] > array[j]) ans = 1;
  else if (array[i] < array[j]) ans = -1;
  return ans;
}

void reverse_int (int *a, int dim) {
  int i=0, j = dim-1, temp;
  while (i < j) {
    temp = a[i];
    a[i] = a[j];
    a[j] = temp;
    i++;
    j--;
  }
}

void check_assignment (int *assigned, int min_k, int max_k, int n) {
  int i, issues = 0, k = max_k - min_k;
  int *clusters = (int *)malloc(sizeof(int)*k);
  fill_array_int(clusters, 0, k);
  for (i = 0; i < n; i++) {
    if ( (assigned[i] < min_k) || (assigned[i] >= max_k) ) {
      // printf("Assignment issue: ");
      // printf("Element %d with invalid assignment %d\n", i, assigned[i]);
      issues++;
    }
    else {
      clusters[min_k + assigned[i]]++;
    }
  }
  printf("%d elements with invalid clusters\n", issues);

  issues = 0;
  for (i = 0; i < k; i++) {
    if (clusters[i] == 0) {
      issues++;
    }
  }
  printf("%d clusters with no elements\n", issues);

  free(clusters);
}

char *concat (char *s1, char *s2) {
  int size = strlen(s1) + strlen(s2) + 1;
  char *ans = (char *)malloc(sizeof(char)*size);
  strcpy(ans, s1);
  strcat(ans, s2);
  ans[size-1] = '\0';
  return ans;
}

double min_value(double *data, int n, int dim) {
  int i;
  double m = data[0];
  for (i = 1; i < n * dim; i++) {
    if (data[i] < m) {
      m = data[i];
    }
  }
  return m;
}

double *get_indices(double *original, int *indices, int n, int dim, int k) {
  int i, idx;
  double *ans = (double *)malloc(sizeof(double) * k * dim);
  for (i = 0; i < k; i++) {
    idx = indices[i];
    copy_to(&original[idx * dim], &ans[i * dim], dim);
  }
  return ans;
}

// Retrieves data (with n rows and fim columns) from a csv file
double* data_from_csv(char* filename, int n, int dim) {
  char ch, *ptr;
  FILE *fp;
  double num;
  int i = 0, j = 0, k = 0;

  char *a = (char *)malloc(sizeof(char) * 30);
  double *data = (double *)malloc(sizeof(double) * n * dim);

  fp = fopen(filename, "r");

  while ( ( ch = getc(fp) ) != EOF) {
    if ( (ch != ',') && (ch != '\n') ) {
      a[i] = ch;
      i++;
    }
    else {
      a[i] = '\0';
      num = atof(a);
      data[j] = num;
      j++;
      i = 0;
    }
  }

  fclose(fp);
  free(a);
  return data;
}

#include "sorting_functions.h"
#include "entropy_functions.h"
#include "printing_functions.h"
#include "random_functions.h"

#endif
