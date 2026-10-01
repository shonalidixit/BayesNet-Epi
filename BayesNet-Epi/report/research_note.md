# BayesNet-Epi — Research Note

## Motivation
Disease transmission in livestock populations can depend strongly on movement connectivity. This project asks how latent network structure, probabilistic risk modelling and partial network observation interact in a controlled synthetic setting.

## Design
A weighted directed network represents movements between synthetic farms. Farms are generated with latent group structure, producing a higher probability of within-group than between-group movement. A stochastic transmission process then propagates infection along directed weighted edges.

## Analysis
Network summaries characterize connectivity and centrality. Greedy modularity provides an exploratory community estimate, compared with the known simulation groups using adjusted Rand index. A Bayesian logistic model relates infection status to standardized weighted in-degree, weighted out-degree and PageRank. Independent Normal(0, 2.5^2) priors are used and posterior uncertainty is approximated locally around the MAP estimate using the inverse Hessian. Finally, repeated induced-node samples quantify distortion of network transitivity.

## Interpretation
Because the ground truth is known, the simulation offers a transparent way to test whether analysis procedures recover designed structure. The aim is methodological demonstration rather than inference about real cattle populations or bovine tuberculosis.

## Limitations and next steps
The network is static and synthetic; transmission lacks recovery and biological heterogeneity; community detection is not a generative SBM; and Laplace inference can underrepresent non-Gaussian posterior uncertainty. A research-grade extension would incorporate temporal movements, spatial information, explicit stochastic block or latent-space models, posterior sampling/variational inference, and—where ethically and legally available—genomic similarity information.
