# make biom
```
biom convert -i genus.xls -o genus.biom --table-type "OTU table" --process-obs-metadata naive --to-hdf5
```
# wilcoxon test
```
wilcox.py -i genus.xls -c LT_S_LT_L -g map.LT_S_LT_L.txt -o LT_S-LT_L.wilcox.genus.xls
```
# Random forest
```
random_forest4key_out_select.pl genus.biom genus.xls map.LT_S_LT_L.txt 0.001
```
# RFCV & ROC
```
# The default value of 'c' is 19,which is changable.
# cross validation
Rscript RFCV.R
## ROC
perl RFCV2ROC.pl genus.xls feature_importance_scores.txt LT_S-LT_L.wilcox.genus.xls map.LT_S_LT_L.txt 16
```
