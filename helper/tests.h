//
// Created by Lucas Murtinho on 31/08/19.
//

#ifndef GITHUBVERSION_TESTS_H

#include "helper_functions.h"

typedef int* (*ClusteringFunc) (double *data, int n, int k, int dim,
                                int m, int max_iter);


void test_single (double *data, int n, int k, int dim, int m,
                      int max_iter, char *dataset, int seed,
                    ClusteringFunc func) {
  FILE *f;
  char *filename = "LC_results.csv";

  if (access(filename, F_OK) == -1 ) {
    f = fopen(filename, "w");
    fprintf(f, "dataset,seed,k,m,time,entropy\n");
  }
  else {
    f = fopen(filename, "a");
  }
  fprintf(f, "%s,%d,%d,%d,", dataset,seed,k,m);
  double *clusters = (double *)malloc(sizeof(double) * k * dim);
  clock_t start, end;
  double t, pe;

  start = clock();
  int *assigned = func(data, n, k, m, dim, max_iter);
  end = clock();
  check_assignment(assigned, 0, k, n);
  maximization(data, assigned, clusters, n, k, dim);
  fprintf(f, "%f,%f\n",
            (double)(end - start) / CLOCKS_PER_SEC,
            partition_entropy(clusters, k, dim));

  free(clusters);
  free(assigned);
  fclose(f);
}

void test_several (double *data, int n, int k, int dim, int m,
                   int max_iter, int times, char *dataset,
                   ClusteringFunc func) {
  int i;
  for (i = 1; i < times+1; i++) {
    srand(i);
    printf("random seed %d:\n", i);
    test_single(data, n, k, dim, m, max_iter, dataset, i, func);
  }
}

void test_full (double *data, int n, int dim, int m, int max_iter,
                    int times, int *k_list, int size_k_list,
                    char *dataset, ClusteringFunc func) {
  int k, i;
  for (i = 0; i < size_k_list; i++) {
    k = k_list[i];
    printf("for %d clusters:\n", k);
    test_several(data, n, k, dim, m, max_iter, times, dataset, func);
  }
}


#define GITHUBVERSION_TESTS_H

#endif //GITHUBVERSION_CLASSES_H
