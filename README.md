# Breeding season management is unlikely to improve population viability of a data-deficient migratory species in decline

### [Kayla L. Davis](https://github.com/davisk93), [Sarah P. Saunders](https://github.com/saund123), Stephanie Beilke, Erin Rowan Ford, Jennifer Fuller, Ava Landgraf, [Elise F. Zipkin](https://zipkinlab.org/)

### 

### Code/Data DOI: TBA upon publication

### Please contact the first author for questions about the code or data: Kayla L. Davis (davisk93@msu.edu).
__________________________________________________________________________________________________________________
## Abstract
A major challenge in conservation is developing effective approaches to mitigate population declines in the face of ongoing environmental change. For migratory species, it is often more feasible to implement management during periods of stationarity, such as the breeding season, when populations are less mobile. However, such management strategies are only successful if the demographic rates they target (e.g., reproductive rates) contribute substantively to population growth. Thus, evaluation of population growth rate sensitivity to variation in demographic parameters is needed to determine the most effective conservation strategies. This is especially true for small and declining populations that require targeted and urgent action to mitigate declines under current and future environmental change. Here, we used a coupled integrated population model-Bayesian population viability analysis (IPM-BPVA) to estimate demographic rates and population viability within the context of climatic and management-related changes for a data-deficient, declining population of black terns in the Upper Midwestern United States. We found that current conservation efforts during the breeding season are unlikely to reverse the declines observed within the last decade (from an average of 307 breeding pairs in 2013 to 54 in 2022). Rather, interventions aimed at increasing adult survival are projected to reduce quasi-extinction probability by 31–48% compared to no additional management or management targeting other rates, depending on the climate scenario.  Our results highlight the importance of enhancing monitoring and management efforts for migratory species during migration and non-breeding periods, which constitute a much larger, and generally riskier, proportion of the annual cycle compared to the breeding season. 

## Data
[SCF-IPM_github.xlsx](https://github.com/davisk93/Davis-et-al_BLTE-IPM-BPVA/blob/main/SCF-IPM_github.xlsx) = Population count data used in the IPM for St. Clair Flats black terns 
[SCF_AHY13to22_update.txt](https://github.com/davisk93/Davis-et-al_BLTE-IPM-BPVA/blob/main/SCF_AHY13to22_update.txt) = Capture histories of black terns marked as adults at St. Clair Flats
[SCF_HY13to22_update.txt](https://github.com/davisk93/Davis-et-al_BLTE-IPM-BPVA/blob/main/SCF_HY13to22_update.txt) = Capture histories of black terns marked as chicks at St. Clair Flats
[YearCovs.csv](https://github.com/davisk93/Davis-et-al_BLTE-IPM-BPVA/blob/main/YearCovs.csv) = Yearly covariates used in the IPM. Only nao_hur (mean value during January-June in year t-1) was supported and used in the final version of the model.
[Hist_NAO_Means.csv](https://github.com/davisk93/Davis-et-al_BLTE-IPM-BPVA/blob/main/Hist_NAO_Means.csv) = Historic values of the North Atlantic Oscillation index during January-June since 1899. We calculated 100-yr (1922-2021) and 10-yr (2012-2021) means of the index to use as priors in the Bayesian population viability analysis (BPVA).

## Code
[BLTE_SCF_Github.R](https://github.com/davisk93/Davis-et-al_BLTE-IPM-BPVA/blob/main/BLTE_SCF_Github.R) = R script containing model code for the IPM and BPVA examining differing management regimes for black terns under two future climate scenarios. 
[Results](https://github.com/davisk93/Davis-et-al_BLTE-IPM-BPVA/tree/main/Results) = This folder contains model results as R data files and code to create figures from the paper. 
