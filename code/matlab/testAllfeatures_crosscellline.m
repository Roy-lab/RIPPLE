% Runs 10fold CV using trained classifiers.
% RIPPLE features.
% Positive examples from Sanyal et al 5C data.

HOME='../../';  % Assume we are running from the MATLAB directory

OUTDIR=sprintf('%s/trainedclassifiers',HOME);
TESTDIR=sprintf('%s/data/5cfeatures',HOME);
train_celltypes={'Helas';'K562';'Gm12878';'H1hesc'}

for c=1:4
		test_celltypes={'Helas';'K562';'Gm12878';'H1hesc'};
		current_train_dir=sprintf('%s/%s/',OUTDIR,train_celltypes{c})
        % load trained forest
		forestfile=sprintf('%s/enhanceronly_forests.mat',current_train_dir);       
		load(forestfile);
        
		for d=1:size(test_celltypes,1)
			if(c==d)
                % don't train/test on same one
				continue
			end
			featurefile=sprintf('%s/%s_enhanceronly_features.txt',TESTDIR,test_celltypes{d})
			current_test_dir=sprintf('%s/../../outputs/%s/',TESTDIR,test_celltypes{d})
            if ~exist(current_test_dir)
                mkdir(current_test_dir);           
            end            
            
            % read testing data
			data=importdata(featurefile);
			fprintf('Test %s with Train %s\n',test_celltypes{d},train_celltypes{c});
			features=data.data;
			rowNames=data.textdata(2:end,1);
			colNames=data.textdata(1,2:end); 
            
			rf_all_actual=features(:,size(features,2));
			folds=length(forests);
            
			labels_soft_all=[];
            sprintf('Feature count:\n');
			fcnt=size(features,2)-1
			featurevals=features(:,1:fcnt);
			randfeaturevals=featurevals(:,randperm(fcnt));
            
            % predict from trained forest
			for t=1:length(forests)
				[final_labels,labels_soft]=predict(forests{t},featurevals);
				labels_soft_all=[labels_soft_all labels_soft(:,2)];	
            end
            
			labels_soft_avg=mean(labels_soft_all,2);
			rf_file=sprintf('%s/rf_train%s_test%s.txt',current_test_dir,train_celltypes{c},test_celltypes{d});
			dlmwrite(rf_file,[labels_soft_avg rf_all_actual],'delimiter','\t')
			rf_withnames_file=sprintf('%s/rfNames_train%s_test%s.txt',current_test_dir,train_celltypes{c},test_celltypes{d});
			fid=fopen(rf_withnames_file,'w');
            sprintf('Verifying sizes of test data labels\n');
			size(rf_all_actual)
			size(labels_soft_avg)
            
            % print out for PR analysis:
            % example, soft label (confidence in positive class), actual
            % class
			for i=1:size(rf_all_actual,1)
				fprintf(fid,'%s\t%d\t%d\n',char(rowNames(i,1)),labels_soft_avg(i,1),rf_all_actual(i,1));
			end
			fclose(fid);
		end
		
end
