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
gdscnet = marginalExp(gdscnet, 0.9, 5, 2);
cclenet = marginalExp(cclenet, 0.9, 5, 2);

%disp('Diffusing mutation over the network...');
%alphadiff = 0.7;
%gdscnet = diffuseMut(gdscnet, net, alphadiff, 1-alphadiff);
%cclenet = diffuseMut(cclenet, net, alphadiff, 1-alphadiff);

maxn = 10;
kfold = 10;

DRUG_NAME = {};
TP = [];
TN = [];
FP = [];
FN = [];
SPEC = [];
PREC = [];
REC = [];

seldrugs = readtable('data/seldrugs.csv', 'Delimiter', ',');
seldrugs.NUMPOS = zeros(height(seldrugs), 1);
seldrugs.NUMNEG = zeros(height(seldrugs), 1);
for i=1:height(seldrugs)
    drug = seldrugs.DRUG_NAME{i};
    thaa = [seldrugs.AATH1(i); seldrugs.AATH2(i)];
    idxd = find(strcmpi(gdscnet.allDrugs, drug));
    idxneg = gdscnet.AAMat(idxd, :) <= thaa(1);
    idxpos = gdscnet.AAMat(idxd, :) >= thaa(2);
    seldrugs.NUMPOS(i) = sum(idxpos);
    seldrugs.NUMNEG(i) = sum(idxneg);
end
conf_th = 0.5;
for i=1:height(seldrugs)
    drug = seldrugs.DRUG_NAME{i};
    idxd = find(strcmpi(gdscnet.allDrugs, drug));
    idxnan = isnan(gdscnet.AAMat(idxd, :));
    thmed = median(gdscnet.AAMat(idxd, ~idxnan));
    idxneg = gdscnet.AAMat(idxd, :) <= thmed;
    idxpos = gdscnet.AAMat(idxd, :) >= thmed;
    fprintf('%s: Neg = %d, Pos = %d\n\n', drug, sum(idxneg), sum(idxpos));
    gdsc_dgnet{i} = drugGeneNetAA(gdscnet, idxd, net, idxpos, idxneg);
    tbl = topGenes(gdsc_dgnet{i}, gdscnet, net, 0.05, 1e-10);
    disp('Doing CV...');

    [pmat, rowlabels, collabels, confmat, predImp, treemdl, predlabel1, predprob1, selmat, sellabels, rf, predlabel2, predprob2] = classifyCV(gdscnet, drug, tbl, maxn, [thmed thmed], kfold);

    writetable(tbl, sprintf('medth/ranking_%s.csv', drug));
    hg = HeatMap(selmat, 'RowLabels', sellabels);
    hg.addTitle(sprintf('Predictor matrix for %s', drug));
    ax = hg.plot;
    colorbar('peer', ax);
    fig = gcf;
    fig.Position(3:4) = 1.5*fig.Position(3:4);
    saveas(gcf, sprintf('medth/%s_predmat.png',drug));
    close all hidden;

    aa = gdscnet.AAMat(idxd, :);
    predlabel = predlabel2;
    predlabel(isnan(aa)) = [];
    predprob = predprob2;
    predprob(isnan(aa)) = [];
    aa(isnan(aa)) = [];
    aap = aa(strcmp(predlabel, 'p') & predprob>=conf_th);
    aan = aa(strcmp(predlabel, 'n') & predprob>=conf_th);
    aau = aa(predprob<conf_th);
    
    
    barx = 0:0.05:1;
    freqp = histc(aap, barx);
    freqn = histc(aan, barx);
    frequ = histc(aau, barx);
    bar(barx, [freqp' freqn' frequ']);
    set(gca, 'XLim', [-0.05 1])
    legend('Positive', 'Negative', 'Undecided');
    title(sprintf('Prediction for %s', drug));
    xlabel('Activity Area');
    ylabel('Frequency');
    print(sprintf('medth/%s_aapred.png', drug), '-dpng');
    close;

    predlabel = predlabel2;
    predprob = predprob2;
    ic50 = gdscnet.IC50Mat(idxd, :);
    predlabel(isnan(ic50)) = [];
    predprob(isnan(ic50)) = [];
    ic50(isnan(ic50)) = [];
    ic50p = ic50(strcmp(predlabel, 'p') & predprob>=conf_th) ;
    ic50n = ic50(strcmp(predlabel, 'n') & predprob>=conf_th);
    ic50u = ic50(predprob<conf_th);
    
   
    mini = min(ic50);
    maxi = max(ic50);
    barx = mini:0.5:maxi;
    freqp = histc(ic50p, barx);
    freqn = histc(ic50n, barx);
    frequ = histc(ic50u, barx);
    bar(barx, [freqp' freqn' frequ']);
    set(gca, 'XLim', [mini-1 maxi+1])
    legend('Positive', 'Negative', 'Undecided');
    title(sprintf('Prediction for %s', drug));
    xlabel('IC50');
    ylabel('Frequency');
    print(sprintf('medth/%s_ic50pred.png', drug), '-dpng');
    close;

    h = treemdl.Impl.viewGraph([], strcmp(treemdl.NodeClass, 'p'),treemdl.PredictorNames, 0,'');
    set(h, 'Position', [0 0 800 600]);
    saveas(h, sprintf('medth/%s_tree.png', drug));
    close(h);
    DRUG_NAME{end+1} = drug;
    TP(end+1) = confmat.tp;
    TN(end+1) = confmat.tn;
    FP(end+1) = confmat.fp;
    FN(end+1) = confmat.fn;
    SPEC(end+1) = confmat.tn/(confmat.tn + confmat.fp);
    PREC(end+1) = confmat.tp/(confmat.tp + confmat.fp);
    REC(end+1) = confmat.tp/(confmat.tp + confmat.fn);
    
    idxdccle = strcmpi(cclenet.allDrugs, drug);
    if(sum(idxdccle)==0)
        continue;
    end
    disp('Predicting on CCLE data...');
    
    predNames = treemdl.PredictorNames;
    predMat = zeros(length(cclenet.cellNames), length(predNames));
    for j=1:length(predNames)
        pp = strsplit(predNames{j}, '-');
        ptype = pp{end};
        genename = strjoin(pp(1:end-1), '-');
        switch(ptype)
            case 'MUT'
                idxgene = strcmpi(cclenet.mutGenes, genename);
                if(sum(idxgene) == 1)
                    predMat(:, j) = cclenet.mutMat(idxgene, :)';
                else
                    warning('Gene %s is not found in the mutation data', genename);
                end
            case 'CNV'
                idxgene = strcmpi(cclenet.cnvGenes, genename);
                if(sum(idxgene) == 1)
                    predMat(:, j) = cclenet.cnvMat(idxgene, :)';
                 else
                    warning('Gene %s is not found in the CNV data', genename);
                end
            case 'GEX'
                idxgene = strcmpi(cclenet.dgexGenes, genename);
                if(sum(idxgene) == 1)
                    predMat(:, j) = cclenet.dgexMat(idxgene, :)';
                else
                    warning('Gene %s is not found in the GEX data', genename);
                end
            otherwise
                error('Error in predictor name %s', predNames{j});
        end
    end
    [predlabel, classprob] = predict(treemdl, predMat);
    
    aa = cclenet.AAMat(idxdccle, :);
    
    predlabel(isnan(aa)) = [];
    classprob(isnan(aa), :) = [];
    aa(isnan(aa)) = [];
    aap = aa(strcmp(predlabel, 'p') & classprob(:, 2)>=conf_th);
    aan = aa(strcmp(predlabel, 'n') & classprob(:, 1)>=conf_th);
    aau = aa(max(classprob, 2)<conf_th);
    
    barx = 0:0.5:8;
    freqp = histc(aap, barx);
    freqn = histc(aan, barx);
    frequ = histc(aau, barx);
    bar(barx, [freqp' freqn' frequ']);
    set(gca, 'XLim', [-1 8])
    legend('Positive', 'Negative', 'Undecided');
    title(sprintf('CCLE Prediction for %s', drug));
    xlabel('Activity Area');
    ylabel('Frequency');
    print(sprintf('medth/%s_predccle.png', drug), '-dpng');
    close all;
    
end

confmattbl = table(DRUG_NAME', TP', TN', FP', FN', SPEC', PREC', REC', 'VariableNames', {'DRUG_NAME'; 'TP';'TN'; 'FP';'FN';'SPEC'; 'PREC'; 'REC'});
writetable(confmattbl, 'medth/medth_confmat.csv');

save -v7.3 medth1 gdsc_dgnet gdscnet net kfold maxn confmattbl; 
