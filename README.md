# net-corr-dim
Code to estimate correlation dimension of complex networks. Developed in MATLAB R2020a.

Associated with the paper\
"Correlation dimension in empirical networks"\
by\
Jack Murdoch Moore, Haiying Wang, Michael Small, Gang Yan, Huijie Yang, and Changgui Gu.\

_Correlation dimension_ D is defined by the power-law
c(s) = s^(D-1),
or the power-law
C(s) = s^D,
where D is correlation dimension, c(s) is correlation (fraction of distinct nodes at distance s), and C(s) is correlation integral (fraction of distinct nodes within network distance s).

_To see examples_ of how to use this code, please run:\
example

_To produce some figures_ from the manuscript, please run:\
show_fit_script\
plot_corr_dim_estimates

_Functions and scripts:_\
compare_corr_dim_est_3.m: Estimate correlation dimension and scaling interval of synthetic networks using different methods and model c(s) = s^(D-1) and save results in folder results-est-dim-3.\
compare_corr_dim_est_4.m: Estimate correlation dimension and scaling interval of synthetic networks using different methods and model C(s) = s^D and save results in folder results-est-dim-4.\
count_distances.m: Return vector of network distances and number of pairs of distinct nodes at each network distance.\
est_corr_dim_3.m: Estimate correlation dimension and scaling interval of a networks using different methods and model c(s) = s^(D-1).\
est_corr_dim_4.m: Estimate correlation dimension and scaling interval of a networks using different methods and model C(s) = s^D.\
example.m: An example to illustrate generation of a synthetic network and estimation of its correlation dimension.\
find_local_minima.m: Find local minima in a vector.\
load_network.m: Load a empirical network from data in folder "networks".\
log_like_3.m: Calculate log-likelihood per observation for model c(s) = s^D.\
log_like_4.m: Calculate log-likelihood per observation for model C(s) = s^D.\
plot_corr_dim_estimates.m: % Plot results of benchmarking, previously saved in folders "results-est-dim-3" and "results-est-dim-4".\
show_fit_func.m: Fit power-law to an interval and illustrate fit and objective function (negative log-likelihood per observation) in figures.\
show_fit_script.m: Illustrate fit and fitting process for empirical or synthetic networks by calling function show_fit_func.m.\
small_world_manhattan.m: Generate lattice* or small world network*.\
small_world_manhattan_lcc.m: Generate lattice* or small world network* and retain only its largest connected component.\

\* Lattices and small world networks are derived from regular $D$-dimensional toroidal lattices defined using a periodic version of the city block (or Manhattan or taxi cab) metric mentioned but not explored in Epidemic dynamics on higher-dimensional small world networks, Applied Mathematics and Computation
421, 126911, by H. Wang, J. M. Moore, M. Small, J. Wang, H. Yang and C. Gu (2022) (associated code at \url{https://github.com/JackMurdochMoore/small-world}).
