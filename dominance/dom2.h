#ifndef DOM2
#define DOM2

#include "../helper/helper_functions.h"
#include "dominance_functions.h"
#include "dom2_functions.h"
#include "dominance_heap.h"
#include "dominance.h"

int* rg2_main(double *data, int n, int k, int dim) {
  int i, j, d;

  // Initial assignment is dominant component (d clusters)
  int *assigned = get_dominant_components(data, n, dim);

  // Ratio is max(a) / sum(a)
  double *ratio = get_ratios(data, n, dim);

  // SORT INDICES OF VECTORS BY ENTROPY
  int *v_indices = largest_elements(ratio, n, n);

  // PARTITION VECTOR INDICES BY FIRST ASSIGNMENT CLUSTER
  int *vinds_per_dom = (int *)malloc(sizeof(int) * dim * n);
  int *vecs_per_dom = (int *)malloc(sizeof(int) * dim);
  first_partition(vinds_per_dom, vecs_per_dom, v_indices, assigned,
    n, dim);

  // Rearrange data according to ratios and dominant components
  int *new_order = order_data(vinds_per_dom, vecs_per_dom, n, dim);
  double *new_data = rearrange(data, new_order, n, dim);

  // Initial assignment of new data: each vector in a cluster
  int *nd_assigned = index_array(n);

  // Get initial entropies of n initial clusters
  double *cluster_ents = get_weighted_entropies(new_data, n, dim);

  int *cutoffs = cum_sum_int(vecs_per_dom, dim);
  int cf_idx = 0;

  int n_clusters = n;
  int h_idx = 0;

  HEAP *heap = (HEAP *)malloc(sizeof(HEAP)*(n-1));
  double *heap_vecs = (double *)malloc(sizeof(double) * (n-1) * dim);
  double *heap_ents = (double *)malloc(sizeof(double)*(n-1));

  i = 0;
  while (i < n-1) {
    if ((cf_idx < dim) && (i == cutoffs[cf_idx]-1)) {
      cf_idx++;
    }
    else {
      copy_to(&new_data[i*dim], &heap_vecs[h_idx*dim], dim);
      add_array(&heap_vecs[h_idx*dim], &new_data[(i+1)*dim], dim);
      heap_ents[i] = weighted_entropy(&heap_vecs[i*dim], dim);

      heap[h_idx] = create_heap_element(new_data, heap_vecs, cluster_ents,
                                          heap_ents, i, n-1, dim);
      h_idx++;
      i++;
    }
  }

  heap_initialize(heap, h_idx);
  HEAP root;
  while (n_clusters > k) {
    root = heap_remove_root(heap, h_idx);
    h_idx--;

    nd_assigned[root.end] = -1;
    add_array(root.ini_vec, root.end_vec, dim);
    cluster_ents[root.ini] = weighted_entropy(root.ini_vec, dim);
    if (root.next >= 0) {
      adjust_heap_element(heap, heap_ents, cluster_ents,
                          root.next, h_idx, dim);
    }
    if (root.prev >= 0) {
      adjust_heap_element(heap, heap_ents, cluster_ents,
                          root.prev, h_idx, dim);
    }
    n_clusters--;
  }

  for (i = 1; i < n; i++) {
    if (nd_assigned[i] == -1) {
      nd_assigned[i] = nd_assigned[i-1];
    }
    else {
      nd_assigned[i] = nd_assigned[i-1] + 1;
    }
  }

  for (i = 0; i < n; i++) {
    assigned[new_order[i]] = nd_assigned[i];
  }

  free(ratio);
  free(v_indices);
  free(vinds_per_dom);
  free(vecs_per_dom);
  free(new_order);
  free(new_data);
  free(nd_assigned);
  free(cutoffs);
  free(heap);
  free(heap_vecs);
  free(cluster_ents);
  free(heap_ents);
  return assigned;
}

int *rg2_clustering(double *data, int n, int k, int dim, double *time) {
  clock_t start = clock();

  int *assigned = (int *)malloc(sizeof(int) * n);
  if (k < dim) {
    double *new_data = transform_data(data, n, k, dim);
    assigned = rg2_main(new_data, n, k, k);
    free(new_data);
  }
  else {
    assigned = rg2_main(data, n, k, dim);
  }

  clock_t end = clock();
  time[0] = (double) (end - start) / CLOCKS_PER_SEC;
  printf("C time for rg2: %f seconds\n", time[0]);
  return assigned;
}

#endif
