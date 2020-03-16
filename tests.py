import ctypes as ct
import numpy as np
import pandas as pd
import scipy.stats as stats
import time

def save_results(df):
    t = time.strftime("%Y_%m_%d_%H_%M")
    df.to_csv("results_"+ t + ".csv")

def update_results(results, data_name, method, k,
                    seed, max_iter, m, t, entropy):
    results['data'].append(data_name)
    results['method'].append(method)
    results['k'].append(k)
    results['seed'].append(seed)
    results['max_iter'].append(max_iter)
    results['m'].append(m)
    results['time'].append(t)
    results['entropy'].append(entropy)

def test_results(data_names, headers, libname, methods,
                    k_list, seed_list, m_list, max_iter_list, verbose=False):
    results = {"data": [],
               "method": [],
               "k": [],
               "seed": [],
               "max_iter": [],
               "m": [],
               "time": [],
               "entropy": []}
    for i in range(len(data_names)):
        data, n, dim = data_from_csv(data_names[i], headers[i])
        if verbose:
            print("testing data set {} with {} rows and {} columns".format(data_names[i], n, dim))
        lib = prepare_lib(libname, methods, n)
        for k in k_list:
            for seed in seed_list:
                if "lc_clustering" in methods:
                    for m in m_list:
                        for max_iter in max_iter_list:
                            if verbose:
                                print("lightcoreset with k = {}, seed = {}, m = {}, max iter = {}".format(k, seed, m, max_iter))
                            pe, t = run_lc_clustering(data, n, dim, k, lib,
                                        m, max_iter, seed)
                            update_results(results, data_names[i],
                                "lc_clustering", k, seed, max_iter, m, t, pe)
                if "rg_clustering" in methods:
                    if verbose:
                        print("ratio-greedy (v1) with k = {}, seed = {}".format(k, seed))
                    pe, t = run_rg_clustering(data, n, dim, k, lib, seed)
                    update_results(results, data_names[i], "rg_clustering", k,
                        seed, np.nan, np.nan, t, pe)
                if "rg2_clustering" in methods:
                    if verbose:
                        print("ratio-greedy with k = {}, seed = {}".format(k, seed))
                    pe, t = run_rg2_clustering(data, n, dim, k, lib, seed)
                    update_results(results, data_names[i], "rg2_clustering", k,
                        seed, np.nan, np.nan, t, pe)
                if "star_clustering" in methods:
                    if verbose:
                        print("star with k = {}, seed = {}".format(k, seed))
                    pe, t = run_star_clustering(data, n, dim, k, lib, seed)
                    update_results(results, data_names[i], "star_clustering", k,
                        seed, np.nan, np.nan, t, pe)
                if "di_clustering" in methods:
                    for max_iter in max_iter_list:
                        if verbose:
                            print("divisive with k = {}, seed = {}, max iter = {}".format(k, seed, max_iter))
                        pe, t = run_di_clustering(data, n, dim, k, lib,
                                    max_iter, seed)
                        update_results(results, data_names[i], "di_clustering",
                            k, seed, max_iter, np.nan, t, pe)
                if "rd_clustering" in methods:
                    if verbose:
                        print("random with k = {}, seed = {}".format(k, seed))
                    pe, t = run_rd_clustering(data, n, dim, k, lib, seed)
                    update_results(results, data_names[i], "rd_clustering", k,
                        seed, np.nan, np.nan, t, pe)
                if "bh_clustering" in methods:
                    for max_iter in max_iter_list:
                        if verbose:
                            print("bregman with k = {}, seed = {}, max iter = {}".format(k, seed, max_iter))
                        pe, t = run_bh_clustering(data, n, dim, k, lib,
                                    max_iter, seed)
                        update_results(results, data_names[i], "bh_clustering",
                            k, seed, max_iter, np.nan, t, pe)
    return results

def get_cluster_vectors(data, assigned):
    return np.array([np.bincount(assigned, weights=data[:,i])
                     for i in range(data.shape[1])]).T

def get_weighted_entropy(v):
    """
    Returns the weighted entropy of a vector.
    Input: 1d array.
    Output: weighted entropy.
    """
    if not v.sum():
        return 0.
    return v.sum() * stats.entropy(v, base=2)

def get_vectors_entropy(cluster_vectors):
    return np.apply_along_axis(get_weighted_entropy, 1, cluster_vectors)

def get_partition_entropy(cluster_vectors):
    return get_vectors_entropy(cluster_vectors).sum()

def analyze_partition(data, n, dim, assigned):
    dm = [i for i in data]
    df = np.reshape(dm, (n, dim))
    cluster_vectors = get_cluster_vectors(df, assigned)
    return get_partition_entropy(cluster_vectors)

def data_from_csv(data_name, header=False):
    """
    Input: root name of a .csv file and whether it has a header.
    Output: the data in the file as a ctypes double array and the original
                dimentions (rows and columns) of the data.
    """
    header = None if (header == False) else "infer"
    file_name = data_name + ".csv"
    df = pd.read_csv(file_name, header=header)
    n, dim = df.shape
    size = df.size
    dm = df.values.reshape((size, 1))
    data = ct.c_double * size
    data = data(*dm)
    return data, n, dim

def prepare_lib(libname, methods, n):
    lib = ct.CDLL(libname)
    for method in methods:
        func = eval("lib.{}".format(method))
        func.restype = np.ctypeslib.ndpointer(dtype=ct.c_int, shape=(n,))
    return lib

def run_lc_clustering(data, n, dim, k, lib,
                        m=3000, max_iter=100, seed=1):
    t = (ct.c_double * 1)()
    assigned = lib.lc_clustering(data, n, k, dim, m, max_iter, t, seed)
    pe = analyze_partition(data, n, dim, assigned)
    # print(pe, t[0])
    return pe, t[0]

def run_rg_clustering(data, n, dim, k, lib, seed=1):
    t = (ct.c_double * 1)()
    assigned = lib.rg_clustering(data, n, k, dim, t)
    pe = analyze_partition(data, n, dim, assigned)
    # print(pe, t[0])
    return pe, t[0]

def run_rg2_clustering(data, n, dim, k, lib, seed=1):
    t = (ct.c_double * 1)()
    assigned = lib.rg2_clustering(data, n, k, dim, t)
    # print(len(np.unique(assigned)))
    pe = analyze_partition(data, n, dim, assigned)
    # print(pe, t[0])
    return pe, t[0]

def run_star_clustering(data, n, dim, k, lib, seed=1):
    t = (ct.c_double * 1)()
    assigned = lib.star_clustering(data, n, k, dim, t)
    # print(len(np.unique(assigned)))
    pe = analyze_partition(data, n, dim, assigned)
    # print(pe, t[0])
    return pe, t[0]

def run_di_clustering(data, n, dim, k, lib,
                        max_iter=100, seed=1):
    t = (ct.c_double * 1)()
    assigned = lib.di_clustering(data, n, k, dim, max_iter, t, seed)
    pe = analyze_partition(data, n, dim, assigned)
    # print(pe, t[0])
    return pe, t[0]

def run_bh_clustering(data, n, dim, k, lib,
                        max_iter=100, seed=1):
    t = (ct.c_double * 1)()
    assigned = lib.bh_clustering(data, n, k, dim, max_iter, t, seed)
    pe = analyze_partition(data, n, dim, assigned)
    # print(pe, t[0])
    return pe, t[0]

def run_rd_clustering(data, n, dim, k, lib, seed=1):
    t = (ct.c_double * 1)()
    assigned = lib.rd_clustering(data, n, k, dim, t, seed)
    pe = analyze_partition(data, n, dim, assigned)
    # print(pe, t[0])
    return pe, t[0]
