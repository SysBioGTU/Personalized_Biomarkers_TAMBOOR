Personalized Prediction of Metabolite Secretion Tendency by the TAMBOOR Algorith
********************************************************************************* 
The TAMBOOR algorithm basically comprises three steps: 
(i) calculation of reaction activity changes between two conditions by mapping transcriptome data on the GEM, 
(ii) calculation of maximum secretion rates of external metabolites by Flux Balance Analysis (FBA), and 
(iii) prediction of number of active reactions required to support the high secretion rates of metabolites, 
ensuring consistency with mapped transcriptome data.  
The codes for each step are given respectively in the files below:
(i) TAMBOOR_step1_mapping.m
(ii) TAMBOOR_step2_max_secretion.m
(iii) TAMBOOR_step3_network_demand.m
Human-GEM (version 1.18.0) with Entrez Gene ID in COBRA version is given in the file "model_ent.mat ".
The Ham’s medium metabolites and corresponding exchange reactions that were used to constraint model in simulations were given the file "rxn.xlsx"

Defining Subgroups of Patients by the Clustering Approach
**********************************************************
The codes for patient grouping with clustering, feature selection, and cluster assignment based on features are given in the file "clustering_patients.R".
Metabolites that can be produced and secreted under given constraints (they are identified in TAMBOOR_step2_max_secretion.m) are given in the file "mets.txt". 
Filtered metabolite lists (dipeptides and tripeptides were filtered out) are given in the file "f_mets.txt". 
