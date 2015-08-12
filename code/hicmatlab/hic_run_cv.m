% Reproduces supplementary experiments from RIPPLE paper.
% Runs 10-fold CV within cell lines GM12878 and K562 from Rao, Huntley HiC.
% Produces 10 output files per cell line - one per fold.

clear all;
K = 10;
cells={'K562','GM12878_combined'}; % Cell lines
enhv=1:1;   % Enhancer stringency mode - only v1 used for RIPPLE supplementary results.

HOME='../../' % assume running from code/hicmatlab directory
RES=sprintf('%s/hic_results', HOME); % result directory
addpath(HOME);

if ~exist(RES)
   mkdir(RES);
end

e=1;    % enhancer type 1 used for all experiments
for c=1:numel(cells)   % cell type
        datfile=sprintf('%s/data/hicfeatures/RaoHuntley_%s_allchr_5kb_enh_v%d_examples.txt', HOME, cells{c}, e);
        if ~exist(datfile)
            continue
        end
     
        
        % output file prefix
        outpref=sprintf('%s/RaoHuntley_%s_allchr_5kb_enh_v%d', RES, cells{c}, e);
        outpref
        
        % skip if this fold done already
        lastout=sprintf('%s_fold_%d.txt',outpref, K);
        if exist(lastout)
            continue
        end

        fulldata = importdata(datfile);
        D=fulldata.data;

        [N, M] = size(D);
        
        N0 = 0; % num pos
        N1 = 0; % num neg
        
        % negatives
        X0 = [];
        X1 = [];
        
        % positives
        Y0 = [];
        Y1 = [];
        
        % split up positives and negatives?
        for i = 1 : N
            if (D(i,end) == 0)
                N0 = N0 + 1;
                X0 = [X0; D(i, 1:end-1)];
                Y0 = [Y0; 0];
            else
                N1 = N1 + 1;
                X1 = [X1; D(i, 1:end-1)];
                Y1 = [Y1; 1];
            end;
        end;

        % bins numbered {1, K}
        % (mod(index, K)+1) will say which bin it does into 

        % for each test fold...
        for i = 1 : K
            Ntrain = 0;
            Ntest = 0;
            clear Xtest;
            clear Ytest;
            clear Xtrain;
            clear Ytrain;
            Xtest = [];
            Ytest = [];
            Xtrain = [];
            Ytrain = [];            
            
            % split up negatives
            for i0 = 1 : N0
                rem=mod(i0, K)+1;       % pick the bin 
                if (rem == i)   % test fold
                    Ytest = [Ytest; 0];
                    Xtest = [Xtest; X0(i0, :)];
                else            % train folds
                    Ytrain = [Ytrain; 0];
                    Xtrain = [Xtrain; X0(i0, :)];
                end;
            end;
            for i1 = 1 : N1
                rem=mod(i1, K)+1;      % pick the bin
                if (rem == i)          % test fold
                    Ytest = [Ytest; 1];
                    Xtest = [Xtest; X1(i1, :)];
                else                   % train folds
                    Ytrain = [Ytrain; 1];
                    Xtrain = [Xtrain; X1(i1, :)];
                end;
            end;
            
            fprintf('Running for test bin %d\n', i);      
            output_file=sprintf('%s_fold_%d.txt',outpref, i);
            run_rf(Xtrain, Ytrain, Xtest, Ytest, output_file, 500, 0.67);            
        end
    end;

