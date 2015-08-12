%Train a Random Forests classifier on each of the 4 cell line data
HOME='../../';  % Assume we are running from the MATLAB directory
DIR=sprintf('%s/data/5cfeatures',HOME);
CELLTYPES={'K562';'Gm12878';'H1hesc';'Helas'};

DISTAL_TYPE='enhanceronly';
for sd=1:1
		SUBDIR=sprintf('%s/',DIR);
for c=1:4
		featurefile=sprintf('%s/%s_%s_features.txt',SUBDIR,CELLTYPES{c},DISTAL_TYPE)
		data=importdata(featurefile);
		features=data.data;
		size(features)
		rowNames=data.textdata(2:end,1);
		colNames=data.textdata(1,2:end); 
		numTrees=500;
		NVarToSample='all';
		FBoot=0.1;
		current_dir=sprintf('%s/outputs/%s/',SUBDIR,CELLTYPES{c})
		mkdir(current_dir);
		currsuff=sprintf('%s/%s_',current_dir,DISTAL_TYPE);
		doRandomForest

 		rf_file=sprintf('%s/rf_%s.txt',current_dir,DISTAL_TYPE)
		dlmwrite(rf_file,[rf_all_labels;rf_all_actual]','delimiter','\t')
		
		rf_withnames_file=sprintf('%s/rfNames_%s.txt',current_dir,DISTAL_TYPE)
		fid=fopen(rf_withnames_file,'w');
		soft_labels=rf_all_labels;
		for i=1:size(rf_all_actual,2)
			fprintf(fid,'%s\t%d\t%d\n',char(rf_all_names(1,i)),soft_labels(1,i),rf_all_actual(1,i));
		end
		fclose(fid)
end
end
