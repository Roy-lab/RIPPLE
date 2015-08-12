The goal of this program is to take a pair of regions (positive and negative pairs) and generate
the features for these pairs. The program is flexible in the sense that it can generate both PRODUCT, CONCAT features
with and without correlation, with and without expression and binary and continuous versions of it.

Usage is:

./genFeatures pairfile location_of_regulatory_genomic_data output_featurefile classlabel[0|1] remove_pairs_withzeros[yes|no] featureencoding[concat|product] RNAseq_expression_file correlation[yes|no] featuretype[binary|continuous]

An example usage is given below:
./genFeatures ../../data/5cregions/H1HESC_enhanceronly_negative.txt exampleinput/H1hesc_examplefeaturefiles.txt exampleoutput/temp_features.txt 0 no concat ../../data/rnaseq/H1hesc_5cpro_exp.txt yes binary
