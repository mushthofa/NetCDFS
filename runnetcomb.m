addpath(genpath('l1ktools-master/'))

disp('Loading Entrez - Gene name mapping...');
entr2name = readtable('data/entrez2name.csv');
entrez2name = containers.Map('keytype', 'double', 'valuetype', 'char');
name2entrez = containers.Map('keytype', 'char', 'valuetype', 'double');
for i=1:height(entr2name)
    entrez2name(entr2name.entrezid(i)) = entr2name.genename{i};
    name2entrez(entr2name.genename{i}) = entr2name.entrezid(i);
end


name2entrez('ABL') = 25;
name2entrez('MLL') = 4297;
name2entrez('EWRS1') = 2130;
name2entrez('FTSJD1') = 55783;
name2entrez('MLL2') = 8085;
name2entrez('MLL3') = 58508;

entrezmap.name2entrez = name2entrez;
entrezmap.entrez2name = entrez2name;

ccle = loadCCLE(entrezmap);
gdsc = loadGDSC();



disp('Loading network...');
% Load network
net = loadNet('net/KEGG-ACSN-HI.csv');

% Restrict data to nodes in the network;

cclenet = restrictNet(ccle, net);
gdscnet = restrictNet(gdsc, net);


idxcommon = cellfun(@(x) ~isempty(find(strcmpi(x, gdsc.allDrugs), 1)), cclenet.allDrugs);
cclenet.allDrugs = cclenet.allDrugs(idxcommon);
cclenet.AAMat = cclenet.AAMat(idxcommon, :);
cclenet.IC50Mat = cclenet.IC50Mat(idxcommon, :);

cclenet.AAMat = cclenet.AAMat/8;

disp('Combining data...');
combnet = combineData(cclenet, gdscnet);

disp('Computing marginals...');
combnet = marginalExp(combnet, 0.9, 5, 2);

clearvars -except entrezmap combnet gdsc ccle net;
save -v7.3 combdata;

tic;
disp('Combined data...');
numall = zeros(length(combnet.allDrugs), 1);
threst = numall;
thresp = numall;
numresp = numall;
numrest = numall;
for i=1:length(combnet.allDrugs)
    disp(combnet.allDrugs{i});
    allprof = combnet.AAMat(i, :);
    allprof  = allprof(~isnan(allprof));
    threst(i) = 1; %quantile(allprof, 0.2);
    thresp(i) = 0; %quantile(allprof, 0.95);
    numrest(i) = sum(allprof<=threst(i));
    numresp(i) = sum(allprof>=thresp(i));
    numall(i) = length(allprof);
    comb_dgnet{i} = drugGeneNet(combnet, combnet.allDrugs{i}, net, 'aa', [threst(i) thresp(i)]);
    save -v7.3 combnet2  comb_dgnet
end
runtime_comb = toc;
save -v7.3 combnet2 comb_dgnet runtime_comb;


