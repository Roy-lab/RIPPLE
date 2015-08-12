% Reproduces supplementary experiments from RIPPLE paper.
% Runs cross-cell-line prediction for HiC GM12878 and K562.
% Train a model on ALL examples for one cell type; test on the other.

clear all;
cells={'K562','GM12878_combined'}; % Cell lines
enhv=1:1;   % Enhancer stringency mode - only v1 used for RIPPLE supplementary results.

HOME='../../' % assume running from code/hicmatlab directory
RES=sprintf('%s/hic_results', HOME); % result directory

if ~exist(RES)
    mkdir(RES);
end

e=1;    % enhancer type 1 only - distal from TSS
clear D1;
clear D2;

% learn forest from each cell type
datfile1=sprintf('%s/data/hicfeatures/RaoHuntley_%s_allchr_5kb_enh_v%d_examples.txt', HOME, cells{1}, e);
datfile2=sprintf('%s/data/hicfeatures/RaoHuntley_%s_allchr_5kb_enh_v%d_examples.txt', HOME, cells{2}, e);

if ~exist(datfile1) || ~exist(datfile2)
    fprintf('Missing data file for enhancer v%d\n', e);
    continue
end

% first cell line
fulldata1 = importdata(datfile1);
D1=fulldata1.data;
X1=D1(:,1:end-1);
Y1=D1(:,end);

% other cell line
fulldata2 = importdata(datfile2);
D2=fulldata2.data;
X2=D2(:,1:end-1);
Y2=D2(:,end);

% assume same features
twohasone=sum(ismember(fulldata1.textdata(1,:),fulldata2.textdata(1,:)));
onehastwo=sum(ismember(fulldata2.textdata(1,:),fulldata1.textdata(1,:)));
if twohasone ~= size(fulldata1.textdata,2) || onehastwo ~=size(fulldata2.textdata,2)
    fprintf('Something weird about features in enhancer v%d\n', e);
    continue
end

fprintf('Cross-cell-line prediction, enhancer version %d\n', e);

% train on 1, test on 2
outtrain1=sprintf('%s/RaoHuntley_predict_%s_from_%s_enh_v%d.list', RES, cells{2}, cells{1}, e);
fprintf('\tPredict %s from %s...\n', cells{2}, cells{1});
if ~exist(outtrain1)
    run_rf(X1,Y1,X2,Y2, outtrain1, 500, 0.67);
end
fprintf('\t%s\n', outtrain1);

% train on 2, test on 1
outtrain2=sprintf('%s/RaoHuntley_predict_%s_from_%s_enh_v%d.list', RES, cells{1}, cells{2}, e);
fprintf('\tPredict %s from %s...\n', cells{1}, cells{2});
if ~exist(outtrain2)
    run_rf(X2,Y2,X1,Y1, outtrain2, 500, 0.67);
end
fprintf('\t%s\n', outtrain2);

