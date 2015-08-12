% Reproduces supplementary experiments from RIPPLE paper.
% Trains RF with RIPPLE features.
% To predict one platform from another: HiC to 5C, and vice versa.

clear all;
K = 10;
cells={'K562','GM12878_combined'}; % Cell lines
enhv=1:1;   % Enhancer type -only 1 used for RIPPLE supplementary results

HOME='../../' % assume running from code/hicmatlab directory
RES=sprintf('%s/hic_results', HOME); % result directory

if ~exist(RES)
    mkdir(RES);
end

for c=1:numel(cells)   % cell type
    for e=enhv  % enhancer type 
        % learn forest from HiC
        datfile=sprintf('%s/data/hicfeatures/RaoHuntley_%s_allchr_5kb_enh_v%d_examples.txt', HOME, cells{c}, e);
               
        if ~exist(datfile)
            continue
        end       
               
        % may need to map between naming conventions
        ocell=cells{c};
        if strcmp(cells{c}, 'GM12878_combined')
            ocell='Gm12878';
        end
        
        % 5C data file
        fivefile=sprintf('%s/data/5cfeatures/%s_enhanceronly_features.txt', HOME, ocell);   
     
        fulldata = importdata(datfile);
        fullCH=fulldata.textdata(1,:);
        D=fulldata.data;
        
        Xhic=D(:,1:end-1);
        Yhic=D(:,end);
               
        % 5c data
        fivedata = importdata(fivefile);
        Dtest=fivedata.data;
        testCH=strrep(fivedata.textdata(1,:), sprintf('%s_', ocell), '');
        testCols=find(ismember(testCH,fullCH));
        % subtract 1 because matrix doesn't have first column
        testCols=testCols-1;        
        Dtest=Dtest(:,testCols);
        
        % split test data into X matrix and Y vector
        X5c=Dtest(:,1:end-1);
        Y5c=Dtest(:,end);        
        
        fprintf('Cell line %s\n', ocell);
             
        % predict 5C from HiC
        output5c=sprintf('%s/%s_predict_5c_from_HiC_enh_v%d.list', RES, ocell, e);
        fprintf('\tPredict 5C from HiC...\n');
        if ~exist(output5c)            
            run_rf(Xhic,Yhic,X5c,Y5c, output5c, 500, 0.67);
        end
        fprintf('\t%s\n', output5c);        
        
         % predict HiC from 5C
        outputHiC=sprintf('%s/%s_predict_HiC_from_5C_enh_v%d.list', RES, ocell, e);
        fprintf('\tPredict HiC from 5C...\n');
        if ~exist(outputHiC)            
            run_rf(X5c,Y5c,Xhic,Yhic, outputHiC, 500, 0.67);
        end
        fprintf('\t%s\n', outputHiC);        
    end
end
