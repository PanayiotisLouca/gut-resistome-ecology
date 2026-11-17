
### Ecological organisation of the human gut resistome

*Panayiotis Louca, Yatong Gao, Joseph Falconer, Mohammadali Khan Mirzaei, Afroditi Kouraki, Jessica Nightingale, Anthony Kelly, Robert Pope, Alessia Visconti, Eleni Kalopedis, Aitor Blanco-Miguez, Davide Bazzani, Mattia Bolzan, Estelle-Emilie Wakeman, Jonas Grove, Ben Auch, Ivan Liachko, Mario Falchi, Alessandra Bandera, Gonzalo Jorquera, Tim D. Spector, Li Deng, Jordana T. Bell, Alexia Karantana, Lindsey A. Edwards, Ana M. Valdes, Cristina Menni*

---


> *This repository contains the analysis code supporting the manuscript:*
> > Louca, P. et al. (2026). *Ecological organisation of the human gut resistome*. [DOI: To be updated]

---

  ## 📂 Repository Structure
```
.
├── 1.Species ARG associations
│   └── 1.1.Species_ARG_mixed_effect_modelling.R
├── 2.Metabolite_ARG_associations
│   ├── 2.1.ARG_metabolites_mixed_effect_modelling.R
│   └── 2.2.ARG_metabolites_mixed_effect_modelling_sensitivity.R
├── 3.Host_phenotype_ARG_associations
│   ├── 3.1.ARG_host_phenotype_mixed_effect_modelling.R
│   └── 3.2. ARG_host_phenotype_mixed_effect_modelling_sensitivity.R
├── 4.Heritability
│   └── 4.1.ARG_heritability.R
└── README.md

```

---

## ℹ️ Repository Information

### 1. Species ARG associations

- `1.1.Species_ARG_mixed_effect_modelling.R`: Tests pairwise associations between 693 prevalent bacterial species and 57 prevalent antimicrobial resistance genes (ARGs) using linear mixed-effects models (CLR-transformed and rank-based inverse-normal transformed abundances). Models are adjusted for age, sex, BMI, sequencing batch, and family relatedness. Used to generate the species–ARG co-abundance network and to identify taxon-specific ecological niches of individual ARGs.

### 2. Metabolite_ARG_associations

- `2.1.ARG_metabolites_mixed_effect_modelling.R`: Tests associations between the 57 prevalent ARGs and selected stool metabolites (secondary bile acids and short-chain fatty acids) using linear mixed-effects models with the same covariate adjustment as above.
- `2.2.ARG_metabolites_mixed_effect_modelling_sensitivity.R`: Sensitivity analysis of the metabolite–ARG associations that additionally adjusts for bacterial Shannon diversity.

### 3. Host_phenotype_ARG_associations

- `3.1.ARG_host_phenotype_mixed_effect_modelling.R`: Tests associations between resistome features (clinically prioritised ARG carriage, ARG richness, ARG Shannon diversity) and host metabolic parameters (triglyceride-glucose index, serum triglycerides, fasting glucose) using linear mixed-effects models.
- `3.2. ARG_host_phenotype_mixed_effect_modelling_sensitivity.R`: Sensitivity analyses that (i) additionally adjust for bacterial Shannon diversity and (ii) model clinically prioritised ARG carriage as a binary outcome using mixed-effects logistic regression.

### 4. Heritability

- `4.1.ARG_heritability.R`: Estimates the contribution of host genetic factors to resistome α-diversity (richness and Shannon diversity) and to the abundance of the 57 prevalent ARGs using twin structural equation models (ACE/AE/E) implemented in the mets package. Models are adjusted for age, sex, BMI, and sequencing batch; the best-fitting model is selected by AIC.

---

## Data Availability

The analyses in this repository were performed using data from the TwinsUK cohort. Access to the underlying datasets is managed by the Department of Twin Research at King’s College London.

Data can be made available to bona fide researchers following the TwinsUK data access procedures:
[twinsuk.ac.uk/resources-for-researchers/access-our-data](https://twinsuk.ac.uk/resources-for-researchers/access-our-data/)

---

## License

This repository is distributed under the MIT License. See the `LICENSE` file for details.