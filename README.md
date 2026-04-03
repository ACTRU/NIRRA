[![DOI](https://zenodo.org/badge/1183506305.svg)](https://doi.org/10.5281/zenodo.19054546)

# NIRRA
Network-Informed Restricted-set Ridge Analysis pipeline for proteomics data

Citation

Tower, W. (2026)

https://doi.org/10.5281/zenodo.19054547

Contact: wtower@mgh.harvard.edu

NIRRA (Network-Informed Restricted-set Ridge Analysis) is a computational framework for discovering biologically meaningful predictive structure in high-dimensional proteomics data. The method leverages prior knowledge of protein–protein interaction networks (as defined by STRING) to construct restricted candidate protein sets that are evaluated using ridge regression models. Rather than treating proteins as independent predictors, NIRRA exploits network topology to identify coordinated protein groups whose combined states are associated with phenotypic outcomes. Predictive sets are subsequently clustered into higher-order modules and decomposed into signed protein interaction motifs, allowing the underlying biological structure of predictive signals to be interpreted. Hierarchical empirical Bayes (HEB) shrinkage is then used to stabilize motif effect estimates across sets with differing representation in the sampled network space. Together, these steps enable robust identification of biologically interpretable proteomic modules and motifs that capture coordinated molecular states associated with complex phenotypes.

NIRRA is designed to recover predictive structure that is often missed by many of the most widely used omics analysis techniques, including differential expression/abundance testing, coexpression-based methods (e.g., WGCNA), and standard machine learning approaches such as random forests and gradient boosting. These approaches typically operate on marginal or correlation-based structure, or capture complex interactions in forms that are difficult to interpret and weakly constrained by biological priors. In contrast, NIRRA explicitly incorporates known biological interaction networks into model construction, enabling the identification of coherent, higher-order feature interactions that are both robustly predictive and mechanistically interpretable. Importantly, rather than committing to a single optimized model as in conventional machine learning, NIRRA systematically evaluates and retains a structured space of network-constrained models, exposing the broader predictive landscape and enabling direct interrogation of the structures that drive predictive performance in interpretable units. As a result, it reveals structured signal that is not readily accessible to conventional pipelines while avoiding the opacity of black-box models.

Pipeline Overview

The NIRRA pipeline consists of the following stages:

1) STRING-based PPI set generation
   - Construction of candidate protein sets using interaction networks derived from STRING.

2) Ridge regression modeling
   - Training predictive models for each protein set against a target phenotype.

3) Set clustering
   - Grouping predictive sets based on protein composition to identify higher-level modules.

4) Module ranking
   - Ranking modules based on predictive performance.

5) Motif-based decomposition and HEB stabilization
   - Decomposition of ridge models into signed protein motifs and stabilization of motif effect estimates using hierarchical empirical Bayes (HEB).

Repository Contents
- Scripts/        Analysis scripts for each pipeline stage
- Best Practices
- Methods Description 
- LICENSE

This project is licensed under the Apache License 2.0. See the LICENSE file for details.

Author

William Tower
