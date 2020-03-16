#ifndef STAR
#define STAR

#include "../helper/helper_functions.h"
#include "../helper/entropy_functions.h"
#include "../dominance/dominance_functions.h"
#include "../dominance/dom2_functions.h"

typedef struct STAR_HEAP_EL {
  double val;
  double *vec;
  int c1;
  int c2;
} STAR_HEAP_EL;

typedef struct STAR_HEAP {
  int n_elements;
  int n;
  int d;
  STAR_HEAP_EL *heap;
  int *heap_els;
  int *n_els;
  int *neighbors;
  int *neighbor_list;
} STAR_HEAP;

void heap_check(STAR_HEAP *heap) {
  int bad = 0;
  int last = heap->n / 2 + 1;
  int s1, s2, i;
  for (i=0; i < last; i++) {
    s1 = i * 2 + 1;
    s2 = s1 + 1;
    if ((s1 < heap->n) && heap->heap[i].val > heap->heap[s1].val) {
      bad++;
    }
    if ((s2 < heap->n) && heap->heap[i].val > heap->heap[s2].val) {
      bad++;
    }
  }
  printf("%d inversions\n", bad);
}

double get_value(double *v1, double *v2, double *v_sum, int d) {
  double ent_sum = weighted_entropy(v_sum, d);
  double ent1 = weighted_entropy(v1, d);
  double ent2 = weighted_entropy(v2, d);
  return ent_sum - ent1 - ent2;
}

void star_heap_add_neighbor(STAR_HEAP *heap, int el, int neighbor) {
  heap->heap_els[el*heap->d + heap->n_els[el]] = heap->n;
  heap->neighbors[el*heap->d + heap->n_els[el]] = neighbor;
  heap->n_els[el]++;
}

void star_heap_add_el(STAR_HEAP *heap, double *data, int i, int j) {
  STAR_HEAP_EL new_el;
  new_el.c1 = i;
  new_el.c2 = j;

  new_el.vec = (double *)malloc(sizeof(double) * heap->d);
  copy_to(&data[i*heap->d], new_el.vec, heap->d);
  add_array(new_el.vec, &data[j*heap->d], heap->d);

  new_el.val = get_value(&data[i*heap->d], &data[j*heap->d],
    new_el.vec, heap->d);

  heap->heap[heap->n] = new_el;

  star_heap_add_neighbor(heap, i, j);
  star_heap_add_neighbor(heap, j, i);

  heap->n++;
}

void star_heap_change_neighbors(STAR_HEAP *heap, int c, int el1, int el2) {
  int i;
  for (i = 0; i < heap->n_els[c]; i++) {
    if (heap->heap_els[c*heap->d + i] == el1) {
      heap->heap_els[c*heap->d + i] = el2;
    }
    else if (heap->heap_els[c*heap->d + i] == el2) {
      heap->heap_els[c*heap->d + i] = el1;
    }
  }
}

void star_heap_swap(STAR_HEAP *heap, int el1, int el2) {
  STAR_HEAP_EL h1 = heap->heap[el1];
  STAR_HEAP_EL h2 = heap->heap[el2];

  star_heap_change_neighbors(heap, h1.c1, el1, el2);
  star_heap_change_neighbors(heap, h1.c2, el1, el2);
  star_heap_change_neighbors(heap, h2.c1, el1, el2);
  star_heap_change_neighbors(heap, h2.c2, el1, el2);

  heap->heap[el1] = h2;
  heap->heap[el2] = h1;
}

int star_heap_select_son(STAR_HEAP *heap, int el) {
  int s1 = 2 * el + 1;
  int s2 = s1 + 1;
  int ans;
  if (s1 >= heap->n) {
    ans = -1;
  }
  else if ((s2 >= heap->n) || (heap->heap[s1].val <= heap->heap[s2].val)) {
    ans = s1;
  }
  else {
    ans = s2;
  }
  return ans;
}

void star_heap_restore_down(STAR_HEAP *heap, int el) {
  int s;
  while (true) {
    s = star_heap_select_son(heap, el);
    if ((s == -1) || (heap->heap[s].val >= heap->heap[el].val)) {
      break;
    }
    else {
      star_heap_swap(heap, el, s);
      el = s;
    }
  }
}

void star_heap_initialize(STAR_HEAP *heap) {
  int i;
  int initial = heap->n / 2;
  for (i = initial; i >= 0; i--) {
    star_heap_restore_down(heap, i);
  }
}

void star_heap_remove_from_neighbors(STAR_HEAP *heap, int i, int to_remove) {
  int n_els = heap->n_els[i];
  int last_el = heap->neighbors[i*heap->d + n_els - 1];

  int j;
  for (j = 0; j < n_els; j++) {
    if (heap->neighbors[i*heap->d + j] == to_remove) {
      heap->neighbors[i*heap->d + j] = last_el;
      heap->neighbors[i*heap->d + n_els-1] = -1;
    }
  }
}

void star_heap_remove_from_heap_els(STAR_HEAP *heap, int i, int to_remove) {
  int n_els = heap->n_els[i];
  int last_el = heap->heap_els[i*heap->d + n_els - 1];

  int j;
  for (j = 0; j < n_els; j++) {
    if (heap->heap_els[i*heap->d + j] == to_remove) {
      heap->heap_els[i*heap->d + j] = last_el;
      heap->heap_els[i*heap->d + n_els-1] = -1;
    }
  }
}

STAR_HEAP_EL star_heap_remove_root(STAR_HEAP *heap) {
  STAR_HEAP_EL root = heap->heap[0];
  star_heap_swap(heap, 0, heap->n - 1);
  heap->n--;
  star_heap_restore_down(heap, 0);

  star_heap_remove_from_neighbors(heap, root.c1, root.c2);
  star_heap_remove_from_neighbors(heap, root.c2, root.c1);
  star_heap_remove_from_heap_els(heap, root.c1, heap->n);
  star_heap_remove_from_heap_els(heap, root.c2, heap->n);

  heap->n_els[root.c1]--;
  heap->n_els[root.c2]--;

  return root;
}

void star_heap_find_new_neighbors(STAR_HEAP *heap, int *new_neighbors,
                                  int *new_n, int keep, int move) {
  fill_array_int(new_neighbors, -1, heap->d);
  new_n[0] = 0;

  heap->neighbor_list[keep] = 1;
  int keep_n = heap->n_els[keep];
  int i, nxt;
  for (i = 0; i < keep_n; i++) {
    nxt = heap->neighbors[keep*heap->d + i];
    if (nxt != move) {
      new_neighbors[new_n[0]] = nxt;
      new_n[0]++;
      heap->neighbor_list[nxt] = 1;
    }
  }

  int move_n = heap->n_els[move];
  for (i = 0; i < move_n; i++) {
    nxt = heap->neighbors[move*heap->d + i];
    if (heap->neighbor_list[nxt] == 0) {
      new_neighbors[new_n[0]] = nxt;
      new_n[0]++;
    }
  }

  heap->neighbor_list[keep] = 0;
  for (i = 0; i < keep_n; i++) {
    nxt = heap->neighbors[keep*heap->d + i];
    heap->neighbor_list[nxt] = 0;
  }
}

void star_heap_restore_up(STAR_HEAP *heap, int el) {
  int p;
  while (el > 0) {
    p = (el - 2 + (el % 2)) / 2;
    if ((p >= 0) && (heap->heap[el].val < heap->heap[p].val)) {
      star_heap_swap(heap, el, p);
      el = p;
    }
    else {
      break;
    }
  }
}

void star_heap_del_el(STAR_HEAP *heap, int el) {
  heap->heap[el].val = -1;
  star_heap_restore_up(heap, el);
  STAR_HEAP_EL root = star_heap_remove_root(heap);
}

void star_heap_remove_neighbors(STAR_HEAP *heap, int el) {
  int i;
  STAR_HEAP_EL root;
  while (heap->n_els[el] > 0) {
    i = heap->heap_els[el*heap->d];
    star_heap_del_el(heap, i);
  }
}

int *star_clustering(double *data, int n, int k, int dim, double *time) {
  clock_t start = clock();

  double *ratios = get_ratios(data, n, dim);
  int *sorted_order = largest_elements(ratios, n, n);
  double *new_data = rearrange(data, sorted_order, n, dim);
  int *dom_comps = get_dominant_components(new_data, n, dim);

  int *clusters = (int *)malloc(sizeof(int) * n);
  fill_array_int(clusters, -1, n);

  int *vinds_per_dom = (int *)malloc(sizeof(int) * dim * n);
  fill_array_int(vinds_per_dom, -1, dim * n);

  int *vecs_per_dom = (int *)malloc(sizeof(int) * dim);
  fill_array_int(vecs_per_dom, 0, dim);

  int i, c;
  for (i = 0; i < n; i++) {
    c = dom_comps[i];
    vinds_per_dom[c*n + vecs_per_dom[c]] = i;
    vecs_per_dom[c]++;
  }

  int d_with_vec = 0;
  for (i = 0; i < dim; i++) {
    if (vecs_per_dom[i] > 0) {
      d_with_vec++;
    }
  }

  int n_heap = n + (d_with_vec * d_with_vec - 3 * d_with_vec) / 2;
  STAR_HEAP *heap = (STAR_HEAP *)malloc(sizeof(STAR_HEAP));
  heap->n = 0;
  heap->n_elements = n;
  heap->d = dim;
  heap->heap = (STAR_HEAP_EL *) malloc(sizeof(STAR_HEAP_EL)*n_heap);

  heap->heap_els = (int *) malloc(sizeof(int) * n * dim);
  fill_array_int(heap->heap_els, -1, n*dim);

  heap->neighbors = (int *) malloc(sizeof(int) * n * dim);
  fill_array_int(heap->neighbors, -1, n * dim);

  heap->n_els = (int *)malloc(sizeof(int) * n);
  fill_array_int(heap->n_els, 0, n);

  heap->neighbor_list = (int *)malloc(sizeof(int) * n);
  fill_array_int(heap->neighbor_list, 0, n);

  int j, m, l, cur, nxt, idx;
  int *nxt_list = (int *)malloc(sizeof(int) * dim);
  int nxt_list_size = 0;

  for (i = 0; i < dim; i++) {
    l = vecs_per_dom[i];
    for (j = 0; j < l; j++) {
      cur = vinds_per_dom[i*n + j];
      if (j < l-1) {
        nxt_list[0] = vinds_per_dom[i*n + j+1];
        nxt_list_size = 1;
      }
      else {
        nxt_list_size = 0;
        for (m = i+1; m < dim; m++) {
          if (vecs_per_dom[m] > 0) {
            idx = vecs_per_dom[m] - 1;
            nxt_list[nxt_list_size] = vinds_per_dom[m*n + idx];
            nxt_list_size++;
          }
        }
      }
      for (m = 0; m < nxt_list_size; m++) {
        nxt = nxt_list[m];
        star_heap_add_el(heap, new_data, cur, nxt);
      }
    }
    nxt_list_size = 0;
  }

  star_heap_initialize(heap);

  c = n;
  int c1, c2, keep, move;
  int *new_n = (int *)malloc(sizeof(int));
  int *new_neighbors = (int *)malloc(sizeof(int) * dim);
  STAR_HEAP_EL root;

  while (c > k) {
    root = star_heap_remove_root(heap);
    c1 = root.c1;
    c2 = root.c2;

    if (c1 > c2) {
      keep = c2;
      move = c1;
    }
    else {
      keep = c1;
      move = c2;
    }
    clusters[move] = keep;

    add_array(&new_data[keep * dim], &new_data[move * dim], dim);
    star_heap_find_new_neighbors(heap, new_neighbors, new_n, keep, move);
    star_heap_remove_neighbors(heap, keep);
    star_heap_remove_neighbors(heap, move);
    for (i = 0; i < new_n[0]; i++) {
      nxt = new_neighbors[i];
      star_heap_add_el(heap, new_data, keep, nxt);
    }
    c--;
  }

  c = 0;
  for (i = 0; i < n; i++) {
    if (clusters[i] == -1) {
      clusters[i] = c;
      c++;
    }
    else {
      clusters[i] = clusters[clusters[i]];
    }
  }

  // CHECK: only needed if new_data used instead of data
  int *ans = (int *)malloc(sizeof(int) * n);

  for (i = 0; i < n; i++) {
    ans[sorted_order[i]] = clusters[i];
  }

  free(ratios);
  free(dom_comps);
  free(sorted_order);
  free(new_data);
  free(clusters);
  free(vinds_per_dom);
  free(vecs_per_dom);

  free(heap->heap);
  free(heap->heap_els);
  free(heap->n_els);
  free(heap->neighbors);
  free(heap->neighbor_list);

  free(nxt_list);
  free(new_n);
  free(new_neighbors);

  clock_t end = clock();
  time[0] = (double) (end - start) / CLOCKS_PER_SEC;
  printf("C time for star: %f seconds\n", time[0]);

  return ans;
}

#endif
