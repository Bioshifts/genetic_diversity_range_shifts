# Genetic diversity impacts climate-induced species range shifts

This repository contains code and data for the study: **"Genetic diversity impacts climate-induced species range shifts"**.

In this study, we test whether species\' genetic diversity---used here as a proxy for adaptive potential---helps explain interspecific variation in the rates at which species are shifting their geographic ranges in response to climate change.

## Background

Anthropogenic climate change is reshaping biodiversity globally, driving species to **either persist in place** (through adaptation or plasticity) or **move in space** (by tracking suitable climates). Species range shifts---especially towards higher latitudes and elevations---have been extensively documented as key responses to climate warming. However, less is known about the **role of evolutionary potential** in shaping these shifts.

Genetic diversity is expected to **enhance adaptive capacity**, influencing both:

-    **Persistence at trailing edges**, by enabling evolutionary rescue under stressful conditions.

-   **Expansion at leading edges**, by facilitating faster adaptation to new environments or enhancing dispersal potential.

While many adaptation studies rely on genomic approaches applied to model species, recent advances in **macrogenetics** enable broader-scale tests of evolutionary hypotheses using large genetic datasets (e.g., Fonseca et al. 2023).

## Research Questions

We ask:

1.   Is genetic diversity associated with the **velocity** of species\' range shifts?

2.  Does this association depend on **climate exposure** (i.e., climate change velocity)?

3.  Does this relationship vary depending on the **position in the range** (leading edge, centroid, trailing edge)?

## Data Sources

We integrated two global datasets:

-    **Range shifts:** 4,673 latitudinal shift estimates across 1,888 species from the [BioShifts database](https://doi.org/10.6084/m9.figshare.7413365.v1), representing climate-induced movements over land and sea.

-   **Genetic diversity:** Intraspecific nucleotide diversity (π) from [Fonseca et al. 2023](https://academic.oup.com/evlett/article/7/5/331/7226644), calculated from mitochondrial and chloroplast DNA via the phylogatR platform.

## Methods Summary

To test whether genetic diversity influences the speed at which species are shifting their ranges due to climate change, we combined two global datasets:

-    **Species range shifts** from the [BioShifts database](https://figshare.com/articles/dataset/BioShifts_a_global_geodatabase_of_climate-induced_species_redistribution_over_land_and_sea/7413365)

-   **Genetic diversity** (nucleotide diversity, π) from [Fonseca et al. 2023](https://academic.oup.com/evlett/article/7/5/331/7226644)

We focused on latitudinal range shifts that occurred in the same direction as climate warming (i.e., poleward shifts).

We used a statistical model (generalized linear mixed-effects model, or GLMM) to test how genetic diversity, climate change velocity, and the position within the species' range (leading edge, centroid, trailing edge) influence range shift rates. We also accounted for potential confounding factors such as differences in study methods, data types, and sampling designs.

To improve model reliability:

-   We used **absolute values** of range shifts and climate change velocities.

-   We **weighted** observations to correct for data imbalance across range positions.

-   We **bootstrapped** the model (12,000 times) to estimate confidence intervals for each effect.

All code used in the study is written in **R 4.4.2**, and is available in this repository.

Key packages:

-    `lme4`, `glmmTMB` -- mixed modeling

-   `MuMIn` -- model fit metrics

-   `bdc`, `rgbif` -- taxonomic harmonization

-   `terra`, `dplyr`, `ggplot2` -- data wrangling and visualization

See the `scripts/` folder for analysis pipelines and the `data/` folder for preprocessed datasets.

