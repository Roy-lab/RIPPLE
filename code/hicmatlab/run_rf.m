function [ success ] = run_rf( Xtrain, Ytrain, Xtest, Ytest, output_file, num_trees, ex_frac )
%run_rf  Trains RF on one data set, tests on another, and prints out the
%results to the output file.
%   Assume Xs have same number of columns.
% ex_frac : fraction of examples used per tree (eg, 0.67)
% num_trees : number of trees (eg, 500)

        numFeat=size(Xtrain,2);
        model = TreeBagger(num_trees, Xtrain, Ytrain, 'FBoot', ex_frac, 'NVarToSample', fix(sqrt(numFeat)));
        [label, label_soft] = predict(model, Xtest);
        clear outputmatrix;
        outputmatrix = [];
        
        % fill in conf from tree, vals from Ytest
        % conf in positive class
        for itest = 1 : size(Ytest,1)
            outputmatrix = [outputmatrix; label_soft(itest, 2), Ytest(itest)]; %, label(itest)];
        end
        dlmwrite(output_file, outputmatrix, 'delimiter', '\t');
        success=1;
end

