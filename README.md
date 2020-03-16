# Code for the Master's thesis "Theoretical and experimental results on information-theoretic clustering"

This repo presents the code needed to reproduce the experiments presented and analyzed in my Master's thesis, "Theoretical and experimental results in information-theoretic clustering". It includes implementations for several algorithms for clustering while minimizing the overall weighted entropy of the resulting partition (the Partition with Minimum Weighted Entropy Problem [PMWEP]), as well as a Jupyter Notebook to reproduce the thesis' results.

The main code of the algorithms was written in C from scratch, to allow for a fair comparison between the algorithms' running times. Python code was used to retrieve/create the data sets used in the experiments, and to run the tests by calling the C code via ctypes.

The following folders and files are included in this repo. Where applicable, I mention the section of the thesis that should be consulted for further information regarding the folder or file.

- `bregman`: folder containing the C implementation of Lloyd's algorithm for clustering[^1] using the ++ initialization[^2] and Kullback-Leibler divergence as the dissimilarity measure between items (section 5.1).
- `coresets`: folder containing the C implementation of the lightweight coreset algorithm[^3] (section 5.3).
- `dhillon`: folder containing the C implementation of the divisive clustering algorithm[^4] (section 5.2).
- `dominance`: folder containing the C implementation of the ratio-greedy algorithm (section 4.3).
- `helper`: folder containing the C functions used to implement two or more of the algorithms analyzed in the thesis.
- `star`: folder containing the C implementation of the Star algorithm (section 4.4).
- `libcluster.c`: C script to create the library of functions to be called with Python.
- `ng20.csv`: the data from the 20 Newsgroups data set used in the experiments (section 6.1.1).
- `poisson.csv`: synthetic data set built in accordance with the description from a previous paper on coreset clustering[^5] (section 6.1.3).
- `Experiments.ipynb`: Jupyter Notebook allowing for the full reproduction of the results presented in the thesis.
- `libclust.so`: the library of C functions retrieved by Python via C types to run the experiments.
- `make_ng20_data.py`: the code needed to generate the `ng20.csv` file.
- `make_poisson_data.py`: the code needed to generate the `poisson.csv` file.
- `make_rcv1_data.py`: the code needed to generate the `rcv1.csv` file.
- `tests.py`: Python function for running tests over the algorithm in the `libclust.so` library.
- `make_libclust.sh`: shell script to generate the `libclust.so` library.

[^1] LLOYD, S. P. "Least squares quantization in PCM." IEEE Transactions on Information Theory, 28 (2):129–137, 1982.

[^2] ARTHUR, D.; VASSILVITSKII, S. "k-means++: The advantages of careful seeding."  In: PROCEEDINGS OF THEEIGHTEENTH ANNUAL ACM-SIAM SYMPOSIUM ON DISCRETE ALGORITHMS, p. 1027–1035, 2007.

[^3] BACHEM, O.; LUCIC, M; KRAUSE, A. "Scalable k-means clustering via lightweight coresets. In: PROCEEDINGS OFTHE 24TH ACM SIGKDD INTERNATIONAL CONFERENCE ON KNOWL-EDGE DISCOVERY & DATA MINING, p. 1119–1127, 2018.

[^4] DHILLON, I.S.; MALLELA, S.; KUMAR, R. "A divisive information-theoretic feature clustering algorithm for text classification. Journal of machine learning research, 3 (Mar):1265–1287, 2003.

[^5] LUCIC, M.; BACHEM, O.; KRAUSE, A. "Strong coresetsfor hard and soft Bregman clustering with applications to exponential family mixtures. In: PROCEEDINGS OF THE 19TH INTERNATIONAL CONFERENCE ON ARTIFICIAL INTELLIGENCE AND STATISTICS, p. 1–9, 2016.
