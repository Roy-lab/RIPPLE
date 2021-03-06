RIPPLE Code Supplement
----------------------
This supplementary material contains code and data files that can be used to
reproduce results.

July 2015
Contact: sroy@biostat.wisc.edu

data/	Regions and features
	Main results based on 5C data from Sanyal et al:
	5cfeatures/	- labeled training examples (feature vectors)
	5cregions/	- positive/negative training regions
	rnaseq/	- cell-line specific expression values (from RNA-seq) for genes matched to promoters

	Supplementary regions based on Rao, Huntley HiC data:
	hicfeatures/ - labeled examples with feature matrices (input to RIPPLE)
	hicregions/ - positive/negative examples
	hicrnaseq/ - expression values from HiC positive/negative regions

code/	
	generatefeatures/   C++ code to generate features (main results)
	
	matlab/ MATLAB code to run/train a classifier using 10 fold cv.
		testAllfeatures_crosscellline.m: load a classifier from trained classifiers and test a feature file to generate predictions. 
		The classifier will output the probability that a given enhancer-promoter pair interact.
		runAllfeatures_crosscellline.m: train a random forests classifier on the 5C data (in data/5cfeatures).
		Results including a trained classifier (a mat file) and the output of training CV will be written in 
		an `outputs' directory where the features were read from. A cell line-specific output directory will be created.
					
	hicmatlab/	MATLAB code to reproduce supplementary assessments of RIPPLE on HiC data for two cell lines (Gm12878 and K562).
		run_rf.m	- function for training/testing RF and printing out PR file
		hic_run_crossplatform.m - predict HiC from 5C, and vice versa
		hic_run_cv.m - 10-fold cross-validation within each cell line for HiC
		hic_run_cross_lines.m - train model on one cell line (HiC) and test on the other

trainedclassifiers/ Trained Random Forests from RIPPLE features (main results)
	Trained Random forests that can be used to predict interactions for new pairs of regions in the format in data/5cfeatures/*enhanceronly.txt feature format

