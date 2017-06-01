addpath(genpath('l1ktools-master/'))
addpath('distributionPlot/');
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
ctrp = loadCTRP(entrezmap);
disp('Loading network...');
% Load network
net = loadNet('net/KEGG-ACSN-HI.csv');

% Restrict data to nodes in the network;

cclenet = restrictNet(ccle, net);
gdscnet = restrictNet(gdsc, net);
ctrpnet = restrictNet(ctrp, net);

disp('Computing marginals...');
gdscnet = marginalExp(gdscnet, 0.9, 5, 2);
cclenet = marginalExp(cclenet, 0.9, 5, 2);
ctrpnet = marginalExp(ctrpnet, 0.9, 5, 2);

mapcg = readtable('data/mapdrug-ctrp-gdsc.txt', 'Delimiter', '\t');
mapcg.CTRP_NAME = upper(mapcg.CTRP_NAME);
for i=1:height(mapcg)
    idxd = strcmp(ctrpnet.allDrugs, mapcg.CTRP_NAME{i});
    ctrpnet.allDrugs{idxd} = mapcg.GDSC_NAME{i};
end

% Copy the data of MLL2, STAG2 and VIM genes from GDSC to CCLE/CTRP between
% the shared cell lines
sharedcell1 = intersect(gdscnet.cellNames, cclenet.cellNames);
idxcGDSC1 = cellfun(@(x) find(strcmpi(gdscnet.cellNames, x)), sharedcell1);
idxcCCLE = cellfun(@(x) find(strcmpi(cclenet.cellNames, x)), sharedcell1);
sharedcell2 = intersect(gdscnet.cellNames, ctrpnet.cellNames);
idxcGDSC2 = cellfun(@(x) find(strcmpi(gdscnet.cellNames, x)), sharedcell2);
idxcCTRP = cellfun(@(x) find(strcmpi(ctrpnet.cellNames, x)), sharedcell2);

idxd = cellfun(@(x) find(ismember(gdscnet.mutGenes, x)), {'MLL2', 'STAG2'});

cclenet.mutGenes(end+1:end+2) = {'MLL2', 'STAG2'};
ctrpnet.mutGenes(end+1:end+2) = {'MLL2', 'STAG2'};

cclenet.mutMat(end+2, end) = 0;
cclenet.mutMat(end-1:end, idxcCCLE) = gdscnet.mutMat(idxd, idxcGDSC1);

ctrpnet.mutMat(end+2, end) = 0;
ctrpnet.mutMat(end-1:end, idxcCTRP) = gdscnet.mutMat(idxd, idxcGDSC2);


% idxd = find(ismember(gdscnet.dgexGenes, {'VIM'}));
% cclenet.dgexGenes(end+1) = {'VIM'};
% cclenet.dgexMat(end+1,:) = 0;
% cclenet.dgexMat(end, idxcCCLE) = gdscnet.dgexMat(idxd, idxcGDSC1);
% ctrpnet.dgexGenes(end+1) = {'VIM'};
% ctrpnet.dgexMat(end+1,:) = 0;
% ctrpnet.dgexMat(end, idxcCTRP) = gdscnet.dgexMat(idxd, idxcGDSC2);


% disp('Diffusing mutation over the network...');
% alphadiff = 0.7;
% gdscnet = diffuseMut(gdscnet, net, alphadiff, 1-alphadiff);
% cclenet = diffuseMut(cclenet, net, alphadiff, 1-alphadiff);

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
seldrugs.NUMCELL = zeros(height(seldrugs), 1);
seldrugs.NUMTEST = zeros(height(seldrugs), 1);
for i=1:height(seldrugs)
    drug = seldrugs.DRUG_NAME{i};
    thaa = [seldrugs.AATH1(i); seldrugs.AATH2(i)];
    idxd = find(strcmpi(gdscnet.allDrugs, drug));
    idxneg = gdscnet.AAMat(idxd, :) <= thaa(1);
    idxpos = gdscnet.AAMat(idxd, :) >= thaa(2);
    seldrugs.NUMPOS(i) = sum(idxpos);
    seldrugs.NUMNEG(i) = sum(idxneg);
    seldrugs.NUMCELL(i) = sum(~isnan(gdscnet.AAMat(idxd, :)));
    seldrugs.NUMTEST(i) = seldrugs.NUMCELL(i) - (seldrugs.NUMPOS(i) + seldrugs.NUMNEG(i));
    
end
conf_th = 0.5;


fpval = fopen('manth/manth_pval.csv', 'w');
ffeats = fopen('manth/allfeats.csv', 'w');

fprintf(fpval, 'DRUG_NAME, NUMCELL, TH1, TH2, NUMPOS, NUMNEG, GDSC_AUC, GDSC_IC50, CCLE_AV, CCLE_AUC, CTRP_AV, CTRP_AUC\n');
fprintf(ffeats, 'DRUG_NAME,FEATURE,SIGN,NETRANK,PREDRANK\n');

nummutp = zeros(height(seldrugs), 1);
nummutn = nummutp;
numcnvp = nummutp;
numcnvn = nummutp;
numgexp = nummutp;
numgexn = nummutp;
numnetp = nummutp;
numnetn = nummutp;
numpredmat = nummutp;

for i=1:height(seldrugs)
    pv_gdsc_auc = nan;
    pv_gdsc_ic50 = nan;
    pv_ccle_auc = nan;
    pv_ctrp_auc = nan;
    
    drug = seldrugs.DRUG_NAME{i};
    thaa = [seldrugs.AATH1(i); seldrugs.AATH2(i)];
    idxd = find(strcmpi(gdscnet.allDrugs, drug));
    idxneg = gdscnet.AAMat(idxd, :) <= thaa(1);
    idxpos = gdscnet.AAMat(idxd, :) >= thaa(2);
    fprintf('%s: Neg = %d, Pos = %d\n\n', drug, sum(idxneg), sum(idxpos));
    if(seldrugs.RUN(i) == 1)
        gdsc_dgnet{i} = drugGeneNetAA(gdscnet, idxd, net, idxpos, idxneg);
    end
    
    tbl = topGenes(gdsc_dgnet{i}, gdscnet, net, 0.05, 1e-10);
    %disp('Doing CV...');
    
    nummutp(i) = sum(cellfun(@(x) ~isempty(x)&~strcmp(x, ''), tbl.MUTP));
    nummutn(i) = sum(cellfun(@(x) ~isempty(x)&~strcmp(x, ''), tbl.MUTN));
    numcnvp(i) = sum(cellfun(@(x) ~isempty(x)&~strcmp(x, ''), tbl.CNVP));
    numcnvn(i) = sum(cellfun(@(x) ~isempty(x)&~strcmp(x, ''), tbl.CNVN));
    numgexp(i) = sum(cellfun(@(x) ~isempty(x)&~strcmp(x, ''), tbl.GEXP));
    numgexn(i) = sum(cellfun(@(x) ~isempty(x)&~strcmp(x, ''), tbl.GEXN));
    numnetp(i) = sum(cellfun(@(x) ~isempty(x)&~strcmp(x, ''), tbl.NETP));
    numnetn(i) = sum(cellfun(@(x) ~isempty(x)&~strcmp(x, ''), tbl.NETN));
   
    [pmat, rowlabels, collabels,  predImp, treemdl, predlabel1, predprob1, selmat, sellabels, rf, predlabel2, predprob2] = classifyCV(gdscnet, drug, tbl, maxn, thaa, kfold);

    
    for j=length(sellabels):-1:1
        fprintf(ffeats, '%s', drug);
        fprintf(ffeats, ',%s', sellabels{j});
        pp = strsplit(sellabels{j}, '-');
        ptype = pp{end};
        genename = strjoin(pp(1:end-1), '-');
        netrankp = find(strcmpi(tbl.NETP, genename));
        if(isempty(netrankp))
            netrankp = inf;
        end
        netrankn = find(strcmpi(tbl.NETN, genename));
        if(isempty(netrankn))
            netrankn = inf;
        end
        switch(ptype)
            case 'MUT'
                drankp = find(strcmpi(tbl.MUTP, genename));
                drankn = find(strcmpi(tbl.MUTN, genename));
            case 'CNV'
                drankp = find(strcmpi(tbl.CNVP, genename));
                drankn = find(strcmpi(tbl.CNVN, genename));
            case 'GEX'
                drankp = find(strcmpi(tbl.GEXP, genename));
                drankn = find(strcmpi(tbl.GEXN, genename));
            otherwise
                error('Error in predictor name %s', predNames{j});
        end
        if(isempty(drankp))
            drankp = inf;
        end
        if(isempty(drankn))
            drankn = inf;
        end
        frankp = min(drankp, netrankp);
        frankn = min(drankn, netrankn); 
        if(frankp < frankn)
            fprintf(ffeats, ',+,%d', frankp);
        else
            fprintf(ffeats, ',-,%d', frankn);
        end
        fprintf(ffeats, ',%d\n', length(sellabels) -j +1);
    end
    
    numpredmat(i) = length(sellabels);
    writetable(tbl, sprintf('manth/ranking_%s.csv', drug));
    hg = HeatMap(selmat, 'RowLabels', sellabels);
    hg.addTitle(sprintf('Predictor matrix for %s', drug));
    ax = hg.plot;
    colorbar('peer', ax);
    fig = gcf;
    fig.Position(3:4) = 1.5*fig.Position(3:4);
    saveas(gcf, sprintf('manth/%s_predmat.png',drug));
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
    
    if(~isempty(aap) && ~isempty(aan))
        pv_gdsc_auc = ranksum(aap, aan);
        
        % Violin plots
        data = [aan'; aap'];
        label = [zeros(length(aan), 1); ones(length(aap), 1)];
        distributionPlot(data, 'groups', label, 'histOpt', 1, 'color', {'y', 'b'})
        text(-0.5, 0.8, sprintf('Rank sum p-val = %s', pv_gdsc_auc));
        xlabel('Predicted response');
        ylabel('Actual Response');
        title(sprintf('%s: AUC Prediction on GDSC', drug));
        print(sprintf('manth/%s_aucvio.png', drug), '-dpng');
    end
    
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
    print(sprintf('manth/%s_aapred.png', drug), '-dpng');
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
    
    if(~isempty(ic50p) && ~isempty(ic50n))
        pv_gdsc_ic50 = ranksum(ic50p, ic50n);
        % Violin plots
        data = [ic50n'; ic50p'];
        label = [zeros(length(ic50n), 1); ones(length(ic50p), 1)];
        distributionPlot(data, 'groups', label, 'histOpt', 1, 'color', {'y', 'b'})
        text(-0.5, 0.9, sprintf('Rank sum p-val = %s', pv_gdsc_ic50));
        xlabel('Predicted response');
        ylabel('Actual Response');
        title(sprintf('%s: IC50 Prediction on GDSC', drug));
        print(sprintf('manth/%s_ic50vio.png', drug), '-dpng');
    end
    
    
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
    print(sprintf('manth/%s_ic50pred.png', drug), '-dpng');
    close;

    h = treemdl.Impl.viewGraph([], strcmp(treemdl.NodeClass, 'p'),treemdl.PredictorNames, 0,'');
    set(h, 'Position', [0 0 800 600]);
    saveas(h, sprintf('manth/%s_tree.png', drug));
    close(h);
%     DRUG_NAME{end+1} = drug;
%     TP(end+1) = confmat.tp;
%     TN(end+1) = confmat.tn;
%     FP(end+1) = confmat.fp;
%     FN(end+1) = confmat.fn;
%     SPEC(end+1) = confmat.tn/(confmat.tn + confmat.fp);
%     PREC(end+1) = confmat.tp/(confmat.tp + confmat.fp);
%     REC(end+1) = confmat.tp/(confmat.tp + confmat.fn);
    
    idxdccle = strcmpi(cclenet.allDrugs, drug);
    allfeatav = 1;
    ccleav = 0;
    if(sum(idxdccle)==1)
        disp('Predicting on CCLE data...');
        predNames = treemdl.PredictorNames;
        predMat = zeros(length(cclenet.cellNames), length(predNames));
        for j=1:length(predNames)
            if(allfeatav == 0)
                break;
            end
            pp = strsplit(predNames{j}, '-');
            ptype = pp{end};
            genename = strjoin(pp(1:end-1), '-');
            switch(ptype)
                case 'MUT'
                    idxgene = strcmpi(cclenet.mutGenes, genename);
                    if(sum(idxgene) == 1)
                        predMat(:, j) = cclenet.mutMat(idxgene, :)';
                    else
                        if(treemdl.predictorImportance(j) > 0)
                            warning('Gene %s is not found in the mutation data', genename);
                            allfeatav = 0;
                        else
                            predMat(:, j) = 0;
                        end
                    end
                case 'CNV'
                    idxgene = strcmpi(cclenet.cnvGenes, genename);
                    if(sum(idxgene) == 1)
                        predMat(:, j) = cclenet.cnvMat(idxgene, :)';
                     else
                        if(treemdl.predictorImportance(j) > 0)
                            warning('Gene %s is not found in the CNV data', genename);
                            allfeatav = 0;
                        else
                            predMat(:, j) = 0;
                        end
                    end
                case 'GEX'
                    idxgene = strcmpi(cclenet.dgexGenes, genename);
                    if(sum(idxgene) == 1)
                        predMat(:, j) = cclenet.dgexMat(idxgene, :)';
                    else
                        if(treemdl.predictorImportance(j) > 0)
                            warning('Gene %s is not found in the GEX data', genename);
                            allfeatav = 0;
                         else
                            predMat(:, j) = 0;
                        end
                    end
                otherwise
                    error('Error in predictor name %s', predNames{j});
            end
        end
        if(allfeatav == 1)
            ccleav = 1;
            [predlabel, classprob] = predict(treemdl, predMat);

            aa = cclenet.AAMat(idxdccle, :);

            predlabel(isnan(aa)) = [];
            classprob(isnan(aa), :) = [];
            aa(isnan(aa)) = [];
            aap = aa(strcmp(predlabel, 'p') & classprob(:, 2)>=conf_th);
            aan = aa(strcmp(predlabel, 'n') & classprob(:, 1)>=conf_th);
            aau = aa(max(classprob, 2)<conf_th);

            if(~isempty(aap) && ~isempty(aan))
                pv_ccle_auc = ranksum(aap, aan);
                % Violin plots
                data = [aan'; aap'];
                label = [zeros(length(aan), 1); ones(length(aap), 1)];
                distributionPlot(data, 'groups', label, 'histOpt', 1, 'color', {'y', 'b'})
                text(-0.5, 0.9, sprintf('Rank sum p-val = %s', pv_ccle_auc));
                xlabel('Predicted response');
                ylabel('Actual Response');
                title(sprintf('%s: AUC Prediction on CCLE', drug));
                print(sprintf('manth/%s_ccleaucvio.png', drug), '-dpng');
            end
    
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
            print(sprintf('manth/%s_predccle.png', drug), '-dpng');
            
        end
        close all;
    end 
    
    
    allfeatav = 1;
    idxdctrp = strcmpi(ctrpnet.allDrugs, drug);
    ctrpav = 0;
    if(sum(idxdctrp)==1)
        disp('Predicting on CTRP data...');
    
        predNames = treemdl.PredictorNames;
        predMat = zeros(length(ctrpnet.cellNames), length(predNames));
        for j=1:length(predNames)
            if(allfeatav == 0)
                break;
            end
            pp = strsplit(predNames{j}, '-');
            ptype = pp{end};
            genename = strjoin(pp(1:end-1), '-');
            switch(ptype)
                case 'MUT'
                    idxgene = strcmpi(ctrpnet.mutGenes, genename);
                    if(sum(idxgene) == 1)
                        predMat(:, j) = ctrpnet.mutMat(idxgene, :)';
                    else
                        if(treemdl.predictorImportance(j)>0)
                            warning('Gene %s is not found in the mutation data', genename);
                            allfeatav = 0;
                        else
                            predMat(:, j) = 0;
                        end
                    end
                case 'CNV'
                    idxgene = strcmpi(ctrpnet.cnvGenes, genename);
                    if(sum(idxgene) == 1)
                        predMat(:, j) = ctrpnet.cnvMat(idxgene, :)';
                     else
                        if(treemdl.predictorImportance(j)>0)
                            warning('Gene %s is not found in the CNV data', genename);
                            allfeatav = 0;
                        else
                            predMat(:, j) = 0;
                        end
                    end
                case 'GEX'
                    idxgene = strcmpi(ctrpnet.dgexGenes, genename);
                    if(sum(idxgene) == 1)
                        predMat(:, j) = ctrpnet.dgexMat(idxgene, :)';
                    else
                        if(treemdl.predictorImportance(j)>0)
                            warning('Gene %s is not found in the GEX data', genename);
                            allfeatav = 0;
                        else
                            predMat(:, j) = 0;
                        end
                    end
                otherwise
                    error('Error in predictor name %s', predNames{j});
            end
        end
        if(allfeatav == 1)
            ctrpav = 1;
            [predlabel, classprob] = predict(treemdl, predMat);

            aa = ctrpnet.AAMat(idxdctrp, :);

            predlabel(isnan(aa)) = [];
            classprob(isnan(aa), :) = [];
            aa(isnan(aa)) = [];
            aap = aa(strcmp(predlabel, 'p') & classprob(:, 2)>=conf_th);
            aan = aa(strcmp(predlabel, 'n') & classprob(:, 1)>=conf_th);
            aau = aa(max(classprob, 2)<conf_th);
            
            if(~isempty(aap) && ~isempty(aan))
                pv_ctrp_auc = ranksum(aap, aan);
                % Violin plots
                data = [aan'; aap'];
                label = [zeros(length(aan), 1); ones(length(aap), 1)];
                distributionPlot(data, 'groups', label, 'histOpt', 1, 'color', {'y', 'b'})
                text(-0.5, 0.6, sprintf('Rank sum p-val = %s', pv_ctrp_auc));
                xlabel('Predicted response');
                ylabel('Actual Response');
                title(sprintf('%s: AUC Prediction on CTRP', drug));
                print(sprintf('manth/%s_ctrpaucvio.png', drug), '-dpng');
            end
    
            barx = 0:0.5:20;
            freqp = histc(aap, barx);
            freqn = histc(aan, barx);
            frequ = histc(aau, barx);
            bar(barx, [freqp' freqn' frequ']);
            set(gca, 'XLim', [-1 21])
            legend('Positive', 'Negative', 'Undecided');
            title(sprintf('CTRP Prediction for %s', drug));
            xlabel('Activity Area');
            ylabel('Frequency');
            print(sprintf('manth/%s_predctrp.png', drug), '-dpng');
            
            close all;
            
        end
       
    end
    fprintf(fpval, '%s, %d, %.2f, %.2f, %d, %d, %f, %f, %d, %f, %d, %f\n', drug, sum(~isnan(gdscnet.AAMat(idxd, :))), thaa(2), thaa(1), sum(idxpos), sum(idxneg), pv_gdsc_auc, pv_gdsc_ic50, ccleav, pv_ccle_auc, ctrpav, pv_ctrp_auc);
end

confmattbl = table(DRUG_NAME', TP', TN', FP', FN', SPEC', PREC', REC', 'VariableNames', {'DRUG_NAME'; 'TP';'TN'; 'FP';'FN';'SPEC'; 'PREC'; 'REC'});
writetable(confmattbl, 'manth/manth_confmat.csv');
fclose(fpval);
fclose(ffeats);


% save -v7.3 manth_final gdsc_dgnet; 
% exit;
