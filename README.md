# WISP

WISP (Weighted In Silico Pathology) is a novel approach to assess intra-tumoral heterogeneity from bulk molecular profiles. Based on predefined pure molecular or histological populations for a particular cancer type, our approach gives a fine description of each tumor in a standardized way. The methodology is based on non-negative least squares regression and quadratic programming optimization for estimating the mixing proportions of distinct populations for a tumoral sample. It can be applied on trancriptomic data or methylation data.

## WISP steps


WISP consists in 3 main steps: 

1- [Calculation of pure population centroid profiles.](#Step1) 

2- [Estimation of mixing proportions of pure populations for each tumoral sample.](#Step2)  

3- [Graphical representations of the mixing proportion estimations for all samples.](#Step3) 


## Citation

For now, please provide a link to this github repository:

<https://github.com/YunaBlum/WISP>

## Install

You may install this package with [devtools]:

[devtools]: https://github.com/hadley/devtools

```{r results='hide', message=FALSE, warning=FALSE}
devtools::install_github("YunaBlum/WISP")

library(WISP)

```


## Illustrative datasets
Illustrative gene expression datasets of presupposed pure histological samples and mixed samples from pleural mesothelioma restricted to specific markers (for computing time reason).

```{r}
data(dataWISP)
```
`dataWISP` is a list containing the following objects:

`datapure`: data frame of gene expression profiles of the presupposed pure histological samples.
```{r}
dataWISP$datapure[1:5,1:5]
##              T036     T025      T022      T004      T016
## ITLN1   11.456054 11.80786 10.205246 12.201666  6.373776
## PRG4    10.643914 13.40374 13.182918 13.249751 11.930280
## PKHD1L1 10.722577 12.25111  7.501498 11.650849  5.084441
## NRG4     9.395325 10.71960 10.827885  9.734774  9.992148
## SLC28A3  6.529778 11.83156  9.828759 11.486970  6.404284
```
`clpure`: vector of histology annotation for the presupposed pure histological samples.
```{r}
head(dataWISP$clpure)
##  T036  T025  T022  T004  T016  T014 
## Epure Epure Epure Epure Epure Epure 
## Levels: Epure normal Spure
```
`data`: data frame of gene expression profiles of mixed histological samples.
```{r}
dataWISP$data[1:5,1:5]
##              T013      T058      T066     T070      T078
## ITLN1    8.178341 11.550232 11.840814 4.877618 11.587651
## PRG4    12.967090 12.756574 13.391280 6.559053 13.054004
## PKHD1L1 11.042867 11.442587 11.789554 4.659080 11.373127
## NRG4     8.818490  8.240092  9.041123 4.707080  8.435392
## SLC28A3 10.638619 10.920845 11.460884 5.814253  9.994315
```
(optional) `histo`: histological annotation for all the samples present in `data` object. 
```{r}
head(dataWISP$histo)
##                       T013                       T058 
## "Epithelioid mesothelioma" "Epithelioid mesothelioma" 
##                       T066                       T070 
## "Epithelioid mesothelioma" "Epithelioid mesothelioma" 
##                       T078                       T080 
## "Epithelioid mesothelioma" "Epithelioid mesothelioma"
```


## Step1: Refine pure samples and calculate pure population centroid profiles
<a name="Step1"></a>


`WISP.getPureCentro` function takes as input a data frame of gene expression restricted to samples considered as pure by the user (in our example `datapure`) and a vector of population names for each sample (in our example `clpure`). If `pureSamples_filtering` is set to TRUE, WISP will perform its sample filtering procedure. More precisely, for each presupposed pure sample, it estimates the proportions of all populations: if the difference between the supposed pure class proportion and the one of the second existing contingent is more that the threshold specified in `pureSamples_deltaTopWeights` argument, the sample is considered as a mixture of populations and is removed. For example, a threshold of 0.6 means that one expects to have at least 60% of the main existing population in a sample to consider it as pure.
If `pureSamples_filtering` is set to FALSE the function will take all the presupposed pure samples to calculate the pure centroid profiles.
For the centroid calculation, the function retrieves the markers that are specific to each population: if you want a particular number of markers per class you can set `nb_markers_selection` to "custom" and `nb_markers_max_perClass` with the desired number of markers. The final selection will be based on ANOVA p-value cutoff (`markers_cutoff_pval_anovatest`) and AUC cutoff (`markers_cutoff_auc`); genes are ranked by their logFC. If `nb_markers_selection` is set to "optim_kappa", the number of markers is chosen so that it optimizes the condition number. 
If you wish particular genes to be present in the final centroid signature you can add a vector of gene names in `add_markers`.



```{r}
resPure = WISP.getPureCentro(data= dataWISP$datapure,cl= dataWISP$clpure, pureSamples_filtering = TRUE, nb_markers_selection = "custom",   nb_markers_max_perClass = 50, markers_cutoff_pval_anovatest = 0.05, markers_pval_anovatest_fdr = TRUE, markers_cutoff_auc = 0.8, pureSamples_deltaTopWeights = 0.6, plot_heatmap = TRUE, add_markers = NULL, col_purePop = c("Epure"="#f46d43","normal"="grey","Spure"="#66c2a5"))
## [1] "Retrieving markers"
## [1] "Centroid calculation"
## [1] "Weight estimation"
## [1] "5 samples were removed T036" "5 samples were removed T016"
## [3] "5 samples were removed T080" "5 samples were removed N018"
## [5] "5 samples were removed N017"
## [1] "Retrieving markers"
## [1] "Centroid calculation"
## [1] "Weight estimation"
## [1] "No sample was removed "
## [1] "Final Centroid calculation (on 17 samples)"

```

The returned object of addition of a heatmap (if `plot_heatmap`=TRUE) is a list with the following objects:   
`genescentro`: centroids of pure populations based on the pure samples.  
`indkept`: pure samples kept by WISP.  
`indremoved`: samples that were not considered as pure and removed by WISP.  
`genesclasses`: list of markers that were used in centroid calculation for each pure population.  

```{r}
lapply(resPure, head)
## $genescentro
##            Epure    Spure    normal
## ITLN1   11.40341 4.847256  7.698763
## PRG4    12.86137 6.867839 10.366323
## PKHD1L1 10.32061 4.643993  7.937093
## NRG4     9.68924 4.802799  5.343581
## SLC28A3 10.69767 5.933601  6.433749
## UPK1B   12.78011 8.357490  8.071400
## 
## $indkept
##  Epure  Epure  Epure  Epure  Epure  Epure 
## "T025" "T022" "T004" "T014" "T131" "T115" 
## 
## $indremoved
##  Epure  Epure  Epure normal normal 
## "T036" "T016" "T080" "N018" "N017" 
## 
## $genesclasses
## $genesclasses$Epure
##  [1] "ITLN1"        "PRG4"         "PKHD1L1"      "NRG4"        
##  [5] "SLC28A3"      "UPK1B"        "LRP2"         "ANXA9"       
##  [9] "HRASLS5"      "TGM1"         "MUC16"        "SOX6"        
## [13] "CALB2"        "PLLP"         "BCO2"         "SGPP2"       
## [17] "FAM153A"      "BNC1"         "CCDC64"       "UPK3B"       
## [21] "LOC100129233" "LRRN4"        "FLRT3"        "CDON"        
## [25] "FBXO22-AS1"   "PEX5L"        "PRR15"        "NPR1"        
## [29] "ALOX12P2"     "MSLN"         "SLC4A4"       "RERG"        
## [33] "ARHGAP44"     "PROCR"        "VTN"          "VIPR2"       
## [37] "LOC100507387" "EFNA5"        "GAL3ST2"      "LOC100507088"
## [41] "WT1"          "RERG-IT1"     "MEIS2"        "KLHL31"      
## [45] "RAET1E"       "RSPO1"        "SH3RF2"       "NMU"         
## [49] "ADH1C"        "SNORA70F"    
## 
## $genesclasses$Spure
##  [1] "ODZ2"         "TNFSF4"       "DKK1"         "EDIL3"       
##  [5] "ANKRD30B"     "RNU5A-8P"     "SERPINE1"     "ADAMTS6"     
##  [9] "ITGA11"       "CDH2"         "ADAMTS5"      "RAB3B"       
## [13] "ANLN"         "CDH13"        "HTR1D"        "MME"         
## [17] "ANKRD1"       "MFAP5"        "CD274"        "FBN2"        
## [21] "MIR548H3"     "RN5S385"      "LOXL2"        "COL5A2"      
## [25] "PDCD1LG2"     "IBSP"         "IGF2BP2"      "TPX2"        
## [29] "MIR3167"      "TNC"          "GREM1"        "CD70"        
## [33] "ENTPD4"       "CXCL9"        "DCDC2"        "ARNTL2"      
## [37] "COL4A1"       "INHBA"        "LOC100507460" "ALPK2"       
## [41] "ADAM12"       "MELK"         "OXTR"         "SLC7A5"      
## [45] "MMP9"         "HS3ST3A1"     "AXL"          "ITGB3"       
## [49] "CENPI"        "ENPP1"       
## 
## $genesclasses$normal
##  [1] "SFTPA1"     "CEACAM6"    "CYP2B7P1"   "SCGB3A2"    "NAPSA"     
##  [6] "SFTPB"      "SFTA3"      "AQP4"       "GKN2"       "ROS1"      
## [11] "WIF1"       "PIGR"       "NDNF"       "ITGB6"      "SLC34A2"   
## [16] "SCN7A"      "CLDN18"     "BMP5"       "CXCL17"     "CES1"      
## [21] "SFTPC"      "LRRK2"      "AGR3"       "C1orf116"   "ATP13A4"   
## [26] "LIPH"       "CACNA2D2"   "CPA3"       "SLC6A14"    "SFTPD"     
## [31] "RTKN2"      "IGLV3-27"   "DMBT1"      "ZNF385B"    "SFTPA2"    
## [36] "MYH11"      "KCNJ15"     "LHFPL3-AS2" "MS4A2"      "ITGA8"     
## [41] "PON3"       "TMPRSS2"    "VEPH1"      "LAMP3"      "HHIP"      
## [46] "IL33"       "TMC5"       "EPCAM"      "KIT"        "MMRN1"
```




## Step2: Estimate mixing proportions of pure populations for each tumoral sample.
<a name="Step2"></a>

`WISP.getWeight` function takes as input a data frame of samples without any prior knowledge and the centroid profiles given by `WISP.getPureClass`. We let the user perform its own normalization method for the expression dataset as well as using unlogged data and centroids which can provide better results. If the centroid result is applied to another technology (example: pure centroid profiles estimated using expression array dataset and then used to estimate weights from RNAseq profiles), we suggest the user to center or scale the data using the `scaling` argument. In order to evaluate the accuracy of the models, the function performs a global F-test for each model and calculates the adjusted R-squared. Based on these two criteria, the function will annotate the sample as "LIMIT" or "OK" in the last column of the output object. Cutoff for the F-test p-value and adjusted R-squared can be set in the arguments `cutoff_gobalFtest` and `Rsquared_cutoff` respectively.


```{r}
resW = WISP.getWeight(dataWISP$data, resPure$genescentro, scaling = c("none", "scale", "center")[1], cutoff_gobalFtest = 0.05, Rsquared_cutoff = 0.2, cutoff_ttest_weights = 0.05)
head(resW)
##      weight.Epure weight.Spure weight.normal dist.Obs.Model Ftest.pvalue
## T013       0.8087       0.1147        0.0766          12.18 1.108085e-50
## T058       0.8475       0.1066        0.0459           8.19 1.024975e-74
## T066       0.9998       0.0000        0.0002           7.33 6.982654e-90
## T070       0.0000       1.0000        0.0000          10.93 1.068890e-60
## T078       0.8389       0.1611        0.0000          11.34 2.113089e-57
## T080       0.6398       0.3070        0.0532          15.95 7.651597e-35
##      Adjusted.R.squared Pvalue.Epure Pvalue.Spure Pvalue.normal
## T013               0.80 2.194218e-55 5.203173e-04   0.003894200
## T058               0.91 5.567254e-80 2.485520e-06   0.009932865
## T066               0.94 1.791295e-96 1.000000e+00   0.989862939
## T070               0.85 1.000000e+00 2.429116e-72   1.000000000
## T078               0.85 4.331416e-61 3.287125e-07   1.000000000
## T080               0.51 3.788155e-32 2.150684e-11   0.121916513
##      weight.Epure.filtered weight.Spure.filtered weight.normal.filtered
## T013                0.8087                0.1147                 0.0766
## T058                0.8475                0.1066                 0.0459
## T066                0.9998                0.0000                 0.0000
## T070                0.0000                1.0000                 0.0000
## T078                0.8389                0.1611                 0.0000
## T080                0.6398                0.3070                 0.0000
##      topWeightedClass deltaTopWeights     CI.Epure     CI.Spure
## T013            Epure          0.6940 [0.75, 0.87] [0.05, 0.18]
## T058            Epure          0.7409  [0.8, 0.89] [0.06, 0.15]
## T066            Epure          0.9996    [0.96, 1]    [0, 0.04]
## T070            Spure          1.0000    [0, 0.06]    [0.94, 1]
## T078            Epure          0.6778  [0.78, 0.9]  [0.1, 0.22]
## T080            Epure          0.3328 [0.56, 0.72] [0.22, 0.39]
##         CI.normal WARNING
## T013 [0.02, 0.13]      OK
## T058 [0.01, 0.08]      OK
## T066    [0, 0.03]      OK
## T070    [0, 0.05]      OK
## T078    [0, 0.05]      OK
## T080    [0, 0.12]      OK
```

<img src="https://github.com/YunaBlum/WISP/raw/master/figures/heatmap_pures.png">

The output is a data frame containing the weight (=proportion) of each pure population, the distance to the model, the F-test p-value, the adjusted R-squared and for each weight the Confidence Interval and the t-test p-value. The last column gives a warning if the model does not meet the criteria described above.

## Step3: Graphical representations of the mixing proportion estimations for all samples.
<a name="Step3"></a>

We propose two different representations. Either a heatmap or a barplot as specified in the argument `plot_type`.
```{r}
WISP.getPlot(resW, annotsup = dataWISP$histo, plot_type = "heatmap", col_purePop = c("Epure"="#f46d43","Spure"="#66c2a5", "normal"="grey"))
```
<img src="https://github.com/YunaBlum/WISP/raw/master/figures/heatmap_weights.png">


The barplot can be divided by the top weight population assignment for each sample or by an optional supplementary annotation. 
```{r}
WISP.getPlot(resW, annotsup = dataWISP$histo, plot_type = "barplot", barplot_split = "by.topweight", col_purePop = c("Epure"="#f46d43","Spure"="#66c2a5", "normal"="grey"))
```
<img src="https://github.com/YunaBlum/WISP/raw/master/figures/barplot_weights.png">

```{r}
WISP.getPlot(resW, annotsup = dataWISP$histo, plot_type = "barplot", barplot_split = "by.annotsup", col_purePop = c("Epure"="#f46d43","Spure"="#66c2a5", "normal"="grey"))
```
<img src="https://github.com/YunaBlum/WISP/raw/master/figures/barplot_weights_order.png">
