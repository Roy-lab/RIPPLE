k=10;
class_counter=1;

rand('seed',1); % We'll make sure to seed before each random action

% Permute the input data
randInd=randperm(size(features,1));
features=features(randInd,:);
rowNames=rowNames(randInd,:);

train_set=features(:,:);

cum_cm=[];

cm=zeros(class_counter);
validation_set_size=floor(size(train_set)/k);

all_labels=[];
all_actual=[];
all_names=[];
all_foldid=[];
forests={};
for fold=1:k
    fold
    if k==1
        train_range=1:floor((1/2)*size(train_set,1))-1;
        validation_range=floor((1/2)*size(train_set,1)):size(train_set,1);
    elseif fold==k % If the number of images isn't evenly divisible by k, tack them on to the last training set
        train_range=1:(fold-1)*validation_set_size;
        validation_range=(fold-1)*validation_set_size+1:size(train_set,1);
    else
        train_range=[1:(fold-1)*validation_set_size,fold*validation_set_size+1:size(train_set,1)];
        validation_range=(fold-1)*validation_set_size+1:fold*validation_set_size;
    end
    
    val_names=rowNames(validation_range,:);
    train_data=train_set(train_range,:);
    validation_data=train_set(validation_range,:);
    max_labels=zeros(length(validation_range),1);
    final_labels=zeros(length(validation_range),1);
    % Handle soft classification
    target_class=1;
    
   %b = TreeBagger(numTrees,train_data(:,1:size(train_data,2)-1),train_data(:,size(train_data,2)),'FBoot',FBoot,'NVarToSample',NVarToSample,'OOBVarImp','on');
    %b = TreeBagger(numTrees,train_data(:,1:size(train_data,2)-1),train_data(:,size(train_data,2)),'FBoot',FBoot,'NVarToSample',NVarToSample,'Prior',[0.3 0.7]);
   b = TreeBagger(numTrees,train_data(:,1:size(train_data,2)-1),train_data(:,size(train_data,2)),'FBoot',FBoot,'NVarToSample',NVarToSample);
    [final_labels,labels_soft]=predict(b,validation_data(:,1:size(validation_data,2)-1));
	labels_soft(1:10,:)
    forests{fold}=b;
    final_labels=str2num(cell2mat(final_labels));
    
   
    actual=validation_data(:,size(validation_data,2));
    these_results=confusionmat(final_labels,actual,'order',[0:1]);

    % Store the results in the cumulative confusion matrix
    [final_labels,actual];
    cm=cm+confusionmat(final_labels,actual,'order',[0:1]);
    
    all_labels=[all_labels,labels_soft(:,2)']; % The posterior for class 1
    all_actual=[all_actual,actual'];
    all_names=[all_names,val_names'];
    all_foldid=[all_foldid;fold*ones(size(validation_range,2),1)];
end

%save(strcat(currsuff,'forests_oob.mat'),'forests','randInd','colNames');
save(strcat(currsuff,'forests.mat'),'forests','randInd','colNames');
%save(strcat(currsuff,'dtree.mat'),'forests','randInd','colNames');

% Accuracy on all the validation data during training
'Validation Results'
cm
accuracy=sum(diag(cm))/sum(sum(cm))

TP=cm(1,1);
FP=cm(1,2);
FN=cm(2,1);
TN=cm(2,2);

TPR=TP/(TP+FN)
recall=TPR
FPR=FP/(TN+FP)
precision=TP/(TP+FP)

%prec_rec(all_labels,all_actual);
rf_all_labels=all_labels;
rf_all_actual=all_actual;
rf_all_names=all_names;
rf_foldid=all_foldid;
