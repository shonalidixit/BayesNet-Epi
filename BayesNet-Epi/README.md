# BayesNet-Epi
### Statistical network modelling and Bayesian risk inference for synthetic livestock disease transmission

An independent research-style project by **Shonali Dixit** demonstrating statistical computing for network epidemiology. It was motivated by methodological questions in livestock movement networks; **all data are synthetic** and this repository is not affiliated with boviSNP or any UCD research project.

## Research questions
1. Can latent community structure be recovered from a directed livestock-movement network?
2. Which network characteristics are associated with simulated infection risk under a probabilistic model?
3. How much structural information is lost when only a fraction of network nodes is sampled?

## Methods
The pipeline generates farms assigned to latent geographic/structural groups and weighted directed animal movements, simulates transmission along movements, calculates network statistics, estimates communities, and fits a Bayesian logistic risk model. The Bayesian analysis uses a Gaussian prior with MAP estimation and a Laplace approximation to quantify posterior uncertainty. A sampling experiment measures how induced node samples distort network transitivity.

## Why these methods?
Network epidemiology is fundamentally relational: transmission opportunity depends not only on farm-level variables but also on movement structure. Community detection probes latent clustering; probabilistic regression quantifies uncertainty in associations between network position and simulated infection; and the sampling experiment asks whether partial observation preserves structural properties.

## Reproduce
```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\\Scripts\\activate
pip install -r requirements.txt
python src/pipeline.py
pytest -q
```

Outputs are written to `data/` and `figures/`.

## Repository structure
- `src/pipeline.py` — end-to-end reproducible simulation and analysis
- `notebooks/` — research walkthrough notebook
- `data/` — generated synthetic outputs
- `figures/` — generated visualisation
- `tests/` — basic reproducibility tests
- `report/` — short research note

## Interpretation and limitations
This is a methodological demonstration, not an epidemiological claim about bovine TB. The simulated transmission mechanism is deliberately simplified; movements are synthetic; community recovery is exploratory rather than a fitted stochastic block model; and the Laplace approximation is not a replacement for a full posterior-sampling workflow. Natural extensions include temporal networks, explicit stochastic block/latent-space likelihoods, MCMC/variational inference, spatial covariates, genomic similarity, and principled network-sampling designs.

## Research integrity
The project intentionally distinguishes *methodological alignment* from *domain results*. It does not use or claim access to boviSNP data, Irish cattle movement data, SNP data, or veterinary records.
