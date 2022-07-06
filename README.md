# ClassifyExosomePathways_ML

The main analyses can be found under Run_Analyses directory. The main analysis scripts are numbered by the order of the analyses and other unnumbered scripts under Run_Analyses and ExtractFeatures are called by the main analysis scripts. Some of the key data are saved under both data and input_feature_data. 

### Sample Preparation and feature extraction

Rscript 01_BinaryExosomeTargets_SamplePreparation.R --outputdir <dir_for_each_binary_classification> --targetfirst <one_exosome_target_category> --targetsecond <another_exosome_target_category>

./02_RunFeatureExtraction.sh -f <dir_for_each_binary_classification>

Rscript 03_Entropyfiltering.R


### Build Random Forest models for binary classifications

./04_RunRandomForest_binaryClass.sh -f <dir_for_each_binary_classification>

Rscript 05_Calc_metrics_importance.R --filedir <outputdir_from_script04>

### Iterative feature selection

./06_Run_IterativeFeatureSelection_zscore0.5.sh -f <dir_for_each_binary_classification> -z <zscore_cutoff>

./07_Calc_metrics_importance_IterativeZscore0.5.sh -f <outputdir_from_script06>

### Build multi-class Random Forest models

Rscript 08_multiClassPrediction_SamplePreparation.R

./09_RunRandomForest_multiClass.sh -f <dir_for_multiclass_classification>

### Run alternative Boruta feature selection

./10_RunRandomForest_Boruta.sh -f <dir_for_each_binary_classification>





