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

disp('Computing marginals...');
cclenet = marginalExp(cclenet, 0.9, 5, 2);
gdscnet = marginalExp(gdscnet, 0.9, 5, 2);


idxcommon = cellfun(@(x) ~isempty(find(strcmpi(x, gdsc.allDrugs), 1)), ccle.allDrugs);
cclenet.allDrugs = cclenet.allDrugs(idxcommon);
cclenet.AAMat = cclenet.AAMat(idxcommon, :)/8;
cclenet.IC50Mat = cclenet.IC50Mat(idxcommon, :);

clearvars -except entrezmap ccle gdsc net cclenet gdscnet;
save -v7.3 dataload;

tic;
disp('CCLE:');
numall = zeros(length(cclenet.allDrugs), 1);
threst = numall;
thresp = numall;
numresp = numall;
numrest = numall;
for i=1:length(cclenet.allDrugs)
    disp(cclenet.allDrugs{i});
    idxpos = cclenet.AAMat(i, :)>=0.7 | cclenet.IC50Mat(i, :)<=0.5;
    idxneg = cclenet.AAMat(i, :)<=0.05 | cclenet.IC50Mat(i, :)>=8;
    ccle_dgnet{i} = drugGeneNet(cclenet, i, net, idxpos, idxneg);
end
runtime_ccle = toc;
save -v7.3 cclenet2 cclenet ccle_dgnet runtime_ccle;


