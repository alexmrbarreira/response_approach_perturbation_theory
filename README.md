# Covariances with the Response Approach

This code computes the covariance matrix of the matter power spectrum using the response approach to perturbation theory. It includes the Gaussian (G) term and the connected non-Gaussian (cNG) term up to 1-loop. It does not include the super-sample covariance (SSC) term.

The response approach was developed in:

- Barreira & Schmidt 2017a, https://arxiv.org/abs/1703.09212
- Barreira & Schmidt 2017b, https://arxiv.org/abs/1705.01092

where it was used to calculate covariance matrices; these are the analysis and ploting scripts for these two papers. 

Consult the above two papers for comparisons against simulations, which are excellent. The response calculations takes however a negligible fraction of the time!

For more applications of the response approach see:

- Barreira, Krause & Schmidt 2018, https://arxiv.org/abs/1711.07467
- Barreira 2019, https://arxiv.org/abs/1901.01243
- Barreira et al. 2019, https://arxiv.org/abs/1904.02070
- Halder & Barreira 2022, https://arxiv.org/abs/2201.05607

### Dependencies

- python (numpy, scipy, matplotlib)

### Code overview

- The files parameters.py and functions.py define global parameters, variables and functions
- The scripts in compute_cov/ execute the covariance calculation
- The scripts in plots/ make plots (figures are stored here too).

### Gallery

Summary of the kinematic regimes in $k_1-k_2$ space

<img src="plots/fig_regimes_v2.png" width="1000" height=auto/>

Stitching of the tree-level standard perturbation theory (SPT) and response results

<img src="plots/fig_tree_transition_mono.png" width="1000" height=auto/>

Contributions along the diagonal

<img src="plots/fig_covariance_diagonal_mono.png" width="600" height=auto/>

Response functions used in the calculation

<img src="plots/fig_responses_eul.png" width="1000" height=auto/>
