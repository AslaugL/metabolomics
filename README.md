# metabolomics 
(Updated) Scripts and pipelines used to analyse metabolomics in my master thesis in biology. Using a different metabolomics dataset found on kaggle (Arnaud Germain, Dinesh K. Barupal, Susan M. Levine, and Maureen R. Hanson, “Comprehensive Metabolomics of ME/CFS.” Kaggle, 2020, doi: 10.34740/KAGGLE/DSV/1542382.)

## Scripts
**tidying** Tidying the sample dataset for analysis<br/>
**univariate** Univariate tests of difference between groups, and script for creating half violin/half boxplots with statistic's from rstatix output. Ex: <br/>
<img src="plots/example_violinboxplot.png" width=50% height=50%><br/>
**count_carbons_doublebonds** Count carbons and double bonds in lipids with naming conventions "class xx:y", "class(xx:y), "class(xx:y\xx:y)" (such as datasets supplied by Zora Bioscience or Biocrates), create heatmap with carbons on y axis and double bonds on x-axis, colored by log2 of fold change between groups <br/>
<img src="plots/LipidCarbonDoublebonds.png" width=50% height=50%><br/>
