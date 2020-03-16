#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 12:31:24 2019

@author: lucasmurtinho
"""

import numpy as np

def poisson_data(n, d, a, b, seed):
    np.random.seed(seed)
    ans = np.zeros((n, d))
    for j in range(d):
        mu = np.random.gamma(a, b)
        ans[:,j] = np.random.poisson(mu, n)
    return ans

N = 10000
D = 50
A = 10
B = 1e3
SEED = 1

dist = poisson_data(N, D, A, B, SEED)

print(dist)