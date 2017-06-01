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


gdsc = loadGDSC();

disp('Loading network...');
% Load network
net = loadNet('net/KEGG-ACSN-HI.csv');

% Restrict data to nodes in the network;


gdscnet = restrictNet(gdsc, net);

disp('Computing marginals...');
gdscnet = marginalExp(gdscnet, 0.9, 5, 2);


clearvars -except entrezmap  gdsc net  gdscnet;

gdsc_ic50th = readtable('data/GDSC_TH.csv');
seldrugs = readtable('data/seldrugs.csv', 'Delimiter', ',');
tic;
disp('GDSC:');
numall = zeros(length(gdscnet.allDrugs), 1);
threst = numall;
thresp = numall;
numresp = numall;
numrest = numall;
for i=1:length(gdscnet.allDrugs)
    drug = gdscnet.allDrugs{i};
   
    if(~sum(strcmpi(seldrugs.DRUG_NAME, drug)))
        continue;
    end
    disp(gdscnet.allDrugs{i});
    idxth = strcmpi(gdsc_ic50th.Drug_Name, gdscnet.allDrugs{i}); 
    if(sum(idxth)==0)
        continue;
    end
    curthic50 = gdsc_ic50th.Threshold(idxth);
    idxpos = gdscnet.AAMat(i, :)>=0.4 | gdscnet.IC50Mat(i, :)<=curthic50;
    idxneg = gdscnet.AAMat(i, :)<=0.1 & gdscnet.IC50Mat(i, :) > curthic50;
%     idxpos = gdscnet.IC50Mat(i, :) <= curthic50;
%     idxneg = gdscnet.IC50Mat(i, :) > curthic50;
    
    gdsc_dgnet{i} = drugGeneNetAA(gdscnet, i, net, idxpos, idxneg);
    save -v7.3 kfgdscnetaa gdscnet gdsc_dgnet;
end
runtime_gdsc = toc;
save -v7.3 kfgdscnetaa gdscnet gdsc_dgnet runtime_gdsc;


