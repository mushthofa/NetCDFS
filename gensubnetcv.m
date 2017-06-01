rng('default');
minleafsize = 10;
maxnumsplit = 50;
eps = 1e-5;
maxiter = 30;
ds = gdscnet;
dstest = cclenet;
allfeats = readtable('manth/allfeats.csv');
cvcclethaa = readtable('cvcclethaa.csv');
fout = fopen(sprintf('subnets/diffkfoldccle_rng_def_%.2f_%d.csv', alpha, maxn), 'w');

accdiff = [];
accnofs = accdiff;
accfs = accdiff;
fprintf(fout, 'DRUG');
%for i=1:length(aacutoffs)
    fprintf(fout, ',tpnofs');
    fprintf(fout, ',fpnofs');
    fprintf(fout, ',tnnofs');
    fprintf(fout, ',fnnofs');
%end
%for i=1:length(aacutoffs)
    fprintf(fout, ',tpfs');
    fprintf(fout, ',fpfs');
    fprintf(fout, ',tnfs');
    fprintf(fout, ',fnfs');
%end
%for i=1:length(aacutoffs)
    fprintf(fout, ',tpdiff');
    fprintf(fout, ',fpdiff');
    fprintf(fout, ',tndiff');
    fprintf(fout, ',fndiff');
%end
fprintf(fout, '\n');

for i=1:height(seldrugs)
    drug = seldrugs.DRUG_NAME{i};
    if(~ismember(drug, cvcclethaa.DRUG))
        continue;
    end
    thaa = [seldrugs.AATH1(i); seldrugs.AATH2(i)];
    idxd = find(strcmpi(gdscnet.allDrugs, drug));
    idxneg = gdscnet.AAMat(idxd, :) <= thaa(1);
    idxpos = gdscnet.AAMat(idxd, :) >= thaa(2);
    
    tbl = readtable(sprintf('manth/ranking_%s.csv', drug));
    posfeats = allfeats.FEATURE(strcmpi(allfeats.DRUG_NAME, drug) & strcmp(allfeats.SIGN, '+'));
    negfeats = allfeats.FEATURE(strcmpi(allfeats.DRUG_NAME, drug) & strcmp(allfeats.SIGN, '-'));
    for j=1:length(posfeats)
        gg = posfeats{j};
        gg = strsplit(gg, '-');
        gg = strjoin(gg(1:end-1), '-');
        posfeats{j} = gg;
    end
    for j=1:length(negfeats)
        gg = negfeats{j};
        gg = strsplit(gg, '-');
        gg = strjoin(gg(1:end-1), '-');
        negfeats{j} = gg;
    end
    mutpgenes = tbl.MUTP;
    mutngenes = tbl.MUTN;
    cnvpgenes = tbl.CNVP;
    if(~iscell(cnvpgenes))
        cnvpgenes = cell(0);
    end
    cnvngenes = tbl.CNVN;
    gexpgenes = tbl.GEXP;
    gexngenes = tbl.GEXN;
    netpgenes = tbl.NETP;
    netngenes = tbl.NETN;


    mutpgenes = trimlist(mutpgenes, maxn);
    mutngenes = trimlist(mutngenes, maxn);
    cnvpgenes = trimlist(cnvpgenes, maxn);
    cnvngenes = trimlist(cnvngenes, maxn);
    gexpgenes = trimlist(gexpgenes, maxn);
    gexngenes = trimlist(gexngenes, maxn);
    netpgenes = trimlist(netpgenes, maxn);
    netngenes = trimlist(netngenes, maxn);

    idxnode = ismember(netpgenes, ds.mutGenes) & ~ismember(netpgenes, mutngenes);
    mutpgenes = union(mutpgenes, netpgenes(idxnode), 'stable');
    idxnode = ismember(netngenes, ds.mutGenes) & ~ismember(netngenes, mutpgenes);
    mutngenes = union(mutngenes, netngenes(idxnode), 'stable');
    idxnode = ismember(netpgenes, ds.cnvGenes) & ~ismember(netpgenes, cnvngenes);
    cnvpgenes = union(cnvpgenes, netpgenes(idxnode), 'stable');
    idxnode = ismember(netngenes, ds.cnvGenes) & ~ismember(netngenes, cnvpgenes);
    cnvngenes = union(cnvngenes, netngenes(idxnode), 'stable');
    idxnode = ismember(netpgenes, ds.dgexGenes) & ~ismember(netpgenes, gexngenes);
    gexpgenes = union(gexpgenes, netpgenes(idxnode), 'stable');
    idxnode = ismember(netngenes, ds.dgexGenes) & ~ismember(netngenes, gexpgenes);
    gexngenes = union(gexngenes, netngenes(idxnode), 'stable');

    allgenes = union(mutpgenes, mutngenes);
    allgenes = union(allgenes, cnvpgenes);
    allgenes = union(allgenes, cnvngenes);
    allgenes = union(allgenes, gexpgenes);
    allgenes = union(allgenes, gexngenes);
    allgenes = union(allgenes, netpgenes);
    allgenes = union(allgenes, netngenes);
    
    nall = length(allgenes);
    
   
    nodes = ismember(net.nodes, allgenes);
    [comp, n, subgraph] = subconncomp(net.mat, nodes);
    uc = unique(comp);
    hh = histc(comp, uc);
    [maxcomp, imax] = max(hh);
    allgeneslab = allgenes;
    for j=1:nall
        gg = allgenes{j};
        idxg = strcmp(ds.mutGenes, gg);
        ssp = sum(sum(ds.mutMat(idxg, idxpos)));
        ssn = sum(sum(ds.mutMat(idxg, idxneg)));
        msg = '';
        if(ssp * ssn > 0)
            msg = strcat(msg, sprintf('(%.0f-%.0f)', ssp/sum(idxpos)*100, ssn/sum(idxneg)*100));
        else
            msg = strcat(msg, '()');
        end
        idxg = strcmp(ds.cnvGenes, gg);
        ssp = sum(sum(ds.cnvMat(idxg, idxpos)~=0));
        ssn = sum(sum(ds.cnvMat(idxg, idxneg)~=0));
        if(ssp * ssn > 0)
            msg = strcat(msg, sprintf('(%.0f-%.0f)', ssp/sum(idxpos)*100, ssn/sum(idxneg)*100));
        else
            msg = strcat(msg, '()');
        end
        allgeneslab{j} = strcat(allgenes{j}, msg);
    end
%     subnet = graph(subgraph, allgeneslab);
%     p = plot(subnet);
%     p.MarkerSize = 10;
%     
%     highlight(p, find(ismember(allgenes, posfeats)), 'NodeColor', 'g')
%     highlight(p, find(ismember(allgenes, negfeats)), 'NodeColor', 'r')
%     
%     fig = gcf;
%     fig.Position(3:4) = 3*fig.Position(3:4);
%     layout(p, 'force');
%     dpfeats = strcat(posfeats, '(+) ');
%     dnfeats = strcat(negfeats, '(-) ');
%     dfeats = strcat(strjoin(dpfeats), ' ', strjoin(dnfeats));
%     msg = sprintf('%s, N = %d, MaxC = %d, %s', upper(drug), nall, maxcomp, dfeats);
%     disp(msg);
%     title(msg);
%     saveas(gcf, sprintf('subnets/%s_subnet.png',drug));
%     close all hidden;
    
    fprintf(fout, '%s', drug);
%     fprintf('%s', drug);
    % Train classification tree without feature selection
    
    featmat = [];
    featmat = [featmat; ds.mutMat(:, [find(idxneg) find(idxpos)])];
    featmat = [featmat; ds.cnvMat(:, [find(idxneg) find(idxpos)])];
    featmat = [featmat; ds.dgexMat(:, [find(idxneg) find(idxpos)])];
    
    featmat = featmat';
    classtarget = [zeros(sum(idxneg), 1); ones(sum(idxpos), 1)];
    pnames = strcat(ds.mutGenes, '-MUT');
    pnames = [pnames; strcat(ds.cnvGenes, '-CNV')];
    pnames = [pnames; strcat(ds.dgexGenes, '-GEX')];

    
%     tt = templateTree('MinLeafSize', minleafsize, 'MaxNumSplit', maxnumsplit);
%     treemdl = fitcensemble(featmat, classtarget, 'CategoricalPredictors', 'all', 'Learners', tt, 'Method', 'GentleBoost', 'NumLearningCycles', ntrees, 'PredictorNames', pnames);
%       
    treemdl = fitctree(featmat, classtarget, 'CategoricalPredictors', 'all', 'PredictorNames', pnames, 'MinLeafSize', minleafsize, 'MaxNumSplit', maxnumsplit);
    usedfeats = pnames(predictorImportance(treemdl)>0);
    try
        pmat = predictorMatrix(dstest, usedfeats);
        featmattest = zeros(length(dstest.cellNames), size(featmat, 2));
        idxuf = cellfun(@(x) find(strcmpi(pnames, x)), usedfeats);
        featmattest(:, idxuf) = pmat';
        classtest = predict(treemdl, featmattest);
        idxdt = strcmpi(dstest.allDrugs, drug);
        aat = dstest.AAMat(idxdt, :)';
        idxdcv = find(strcmpi(cvcclethaa.DRUG, drug));
        thaas0 = cvcclethaa.THLO(idxdcv);
        thaas1 = cvcclethaa.THHI(idxdcv);
%         for j=1:length(thaas0)
            neglabel = aat <= thaas0;
            poslabel = aat >= thaas1;
            tp = sum(classtest == 1 & poslabel == 1);
            fp = sum(classtest == 1 & neglabel == 1);
            tn = sum(classtest == 0 & neglabel == 1);
            fn = sum(classtest == 0 & poslabel == 1);
            fprintf(fout, ',%d,%d,%d,%d', tp, fp, tn, fn);
%         end
        accnofs(end+1) = (tp+tn)/(tp+fp+tn+fn);
    catch err
        %fprintf('%s: %s', drug, err.message);
%         for j=1:length(aacutoffs)
            fprintf(fout, ',NaN,NaN,NaN,NaN');
            accnofs(end+1) = nan;
%         end
    end
    
    
    
    % Train classification tree WITH feature selection
    
    idxmut = ismember(ds.mutGenes, allgenes);
    idxcnv = ismember(ds.cnvGenes, allgenes);
    idxgex = ismember(ds.dgexGenes, allgenes);
    
    featmat = [];
    featmat = [featmat; ds.mutMat(idxmut, [find(idxneg) find(idxpos)])];
    featmat = [featmat; ds.cnvMat(idxcnv, [find(idxneg) find(idxpos)])];
    featmat = [featmat; ds.dgexMat(idxgex, [find(idxneg) find(idxpos)])];
    
  
    featmat = featmat';
    classtarget = [zeros(sum(idxneg), 1); ones(sum(idxpos), 1)];
    pnames = strcat(ds.mutGenes(idxmut), '-MUT');
    pnames = [pnames; strcat(ds.cnvGenes(idxcnv), '-CNV')];
    pnames = [pnames; strcat(ds.dgexGenes(idxgex), '-GEX')];
    
%     tt = templateTree('MinLeafSize', minleafsize, 'MaxNumSplit', maxnumsplit);
%     treemdl = fitcensemble(featmat, classtarget, 'CategoricalPredictors', 'all', 'Learners', tt, 'Method', 'GentleBoost', 'NumLearningCycles', ntrees, 'PredictorNames', pnames);
      
    treemdl = fitctree(featmat, classtarget, 'CategoricalPredictors', 'all', 'PredictorNames', pnames, 'MinLeafSize', minleafsize, 'MaxNumSplit', maxnumsplit);
    usedfeats = pnames(predictorImportance(treemdl)>0);
    try
        pmat = predictorMatrix(dstest, usedfeats);
        featmattest = zeros(length(dstest.cellNames), size(featmat, 2));
        idxuf = cellfun(@(x) find(strcmpi(pnames, x)), usedfeats);
        featmattest(:, idxuf) = pmat';
        classtest = predict(treemdl, featmattest);
        idxdt = strcmpi(dstest.allDrugs, drug);
        aat = dstest.AAMat(idxdt, :)';
        idxdcv = find(strcmpi(cvcclethaa.DRUG, drug));
        thaas0 = cvcclethaa.THLO(idxdcv);
        thaas1 = cvcclethaa.THHI(idxdcv);
%         for j=1:length(thaas0)
            neglabel = aat <= thaas0;
            poslabel = aat >= thaas1;
            tp = sum(classtest == 1 & poslabel == 1);
            fp = sum(classtest == 1 & neglabel == 1);
            tn = sum(classtest == 0 & neglabel == 1);
            fn = sum(classtest == 0 & poslabel == 1);
            fprintf(fout, ',%d,%d,%d,%d', tp, fp, tn, fn);
%         end
            accfs(end+1) = (tp+tn)/(tp+fp+tn+fn);
    catch err
        %fprintf('%s: %s', drug, err.message);
%         for j=1:length(aacutoffs)
            fprintf(fout, ',NaN,NaN,NaN,NaN');
            accfs(end+1) = nan;
%         end
    end
    
    
    
    % Classification tree but diffuse mutation + CNV first
    
    mutgenes = union(mutpgenes, mutngenes);
    idxmut = cellfun(@(x) find(ismember(ds.mutGenes, x)), mutgenes);
    mutmat = ds.mutMat(idxmut, [find(idxneg) find(idxpos)]);
    
    
    idxmutnet = cell2mat(values(net.name2node, mutgenes));
    adjmat = net.mat(idxmutnet, idxmutnet);
    adjmat = single((adjmat + adjmat')>0);
    
    
    nSample = sum(idxpos)+sum(idxneg);
    diffmutMat = netSmooth(mutmat', adjmat, alpha, eps, maxiter)';   
    diffmutMat = double(diffmutMat >= (1-alpha));
    
    % diffuse amplifications
    cnvgenes = union(cnvpgenes, cnvngenes);
    idxcnv = cellfun(@(x) find(ismember(ds.cnvGenes, x)), cnvgenes);
    cnvmat = single(ds.cnvMat(idxcnv, [find(idxneg) find(idxpos)]) > 0); % Find only 1's
    
    idxcnvnet = cell2mat(values(net.name2node, cnvgenes));
    adjmat = net.mat(idxcnvnet, idxcnvnet);
    adjmat = single((adjmat + adjmat')>0);
    
    diffAmpMat = netSmooth(cnvmat', adjmat, alpha, eps, maxiter)';

    % Diffuse Deletions
    cnvgenes = union(cnvpgenes, cnvngenes);
    idxcnv = cellfun(@(x) find(ismember(ds.cnvGenes, x)), cnvgenes);
    cnvmat = single(ds.cnvMat(idxcnv, [find(idxneg) find(idxpos)]) < 0); % Find only -1's
    
    idxcnvnet = cell2mat(values(net.name2node, cnvgenes));
    adjmat = net.mat(idxcnvnet, idxcnvnet);
    adjmat = single((adjmat + adjmat')>0);
    
    diffDelMat = netSmooth(cnvmat', adjmat, alpha, eps, maxiter)';
    
    diffCNVMat = double(abs(diffAmpMat - diffDelMat) >= (1-alpha)) .* sign((diffAmpMat - diffDelMat));
    
    featmat = [];
    featmat = [featmat; diffmutMat]; 
    featmat = [featmat; diffCNVMat];
    featmat = [featmat; ds.dgexMat(idxgex, [find(idxneg) find(idxpos)])];
    featmat = featmat';
   
    pnames = strcat(mutgenes, '-MUT');
    pnames = [pnames; strcat(ds.cnvGenes(idxcnv), '-CNV')];
    pnames = [pnames; strcat(ds.dgexGenes(idxgex), '-GEX')];

%     tt = templateTree('MinLeafSize', minleafsize, 'MaxNumSplit', maxnumsplit);
%     treemdl = fitcensemble(featmat, classtarget, 'CategoricalPredictors', 'all', 'Learners', tt, 'Method', 'GentleBoost', 'NumLearningCycles', ntrees, 'PredictorNames', pnames);
%        
    treemdl = fitctree(featmat, classtarget, 'CategoricalPredictors', 'all', 'PredictorNames', pnames, 'MinLeafSize', minleafsize, 'MaxNumSplit', maxnumsplit);
    usedfeats = pnames(predictorImportance(treemdl)>0);
    try

        % Get and Diffuse mutations
        idxmut = cellfun(@(x) find(ismember(dstest.mutGenes, x)), mutgenes, 'UniformOutput', 0);
        idxnav = cellfun(@(x) isempty(x), idxmut);
        idxmutav = cell2mat(idxmut(~idxnav));
        mutmattest = dstest.mutMat(idxmutav, :);


        idxmutnet = cell2mat(values(net.name2node, mutgenes(~idxnav)));
        adjmat = net.mat(idxmutnet, idxmutnet);
        adjmat = single((adjmat + adjmat')>0);

        diffmutMat = netSmooth(mutmattest', adjmat, alpha, eps, maxiter)';   
        diffmutMat = double(diffmutMat >= (1-alpha));
        
        
        featmutmat = zeros(length(mutgenes), length(dstest.cellNames));
        featmutmat(~idxnav, :) = diffmutMat;
        
        
        % Get CNV data
        cnvgenes = ds.cnvGenes(idxcnv);
        idxcnv = cellfun(@(x) find(ismember(dstest.cnvGenes, x)), cnvgenes, 'UniformOutput', 0);
        idxnav = cellfun(@(x) isempty(x), idxcnv);
        idxcnvav = cell2mat(idxcnv(~idxnav));
        cnvmattest = dstest.cnvMat(idxcnvav, :);
        
        idxcnvnet = cell2mat(values(net.name2node, cnvgenes(~idxnav)));
        adjmat = net.mat(idxcnvnet, idxcnvnet);
        adjmat = single((adjmat + adjmat')>0);
        
        cnvmattestamp = double(cnvmattest>0);      % Find amplifications
        diffAmpMat = netSmooth(cnvmattestamp', adjmat, alpha, eps, maxiter)';   
        cnvmattestdel = double(cnvmattest<0);      % Find deletions
        diffDelMat = netSmooth(cnvmattestdel', adjmat, alpha, eps, maxiter)';   
        
        diffCNVMat = double(abs(diffAmpMat - diffDelMat) >= (1-alpha)) .* sign((diffAmpMat - diffDelMat));
        
        featcnvmat = zeros(length(cnvgenes), length(dstest.cellNames));
        featcnvmat(~idxnav, :) = diffCNVMat;
        
        % Get GEX data
        gexgenes = ds.dgexGenes(idxgex);
        idxgex = cellfun(@(x) find(ismember(dstest.dgexGenes, x)), gexgenes, 'UniformOutput', 0);
        idxnav = cellfun(@(x) isempty(x), idxgex);
        idxgexav = cell2mat(idxgex(~idxnav));
        gexmattest = dstest.dgexMat(idxgexav, :);
        featgexmat = zeros(length(gexgenes), length(dstest.cellNames));
        featgexmat(~idxnav, :) = gexmattest;

        
        featmattest = [];
        featmattest = [featmattest; featmutmat]; 
        featmattest = [featmattest; featcnvmat];
        featmattest = [featmattest; featgexmat];
        featmattest = featmattest';
    
    
        classtest = predict(treemdl, featmattest);
        idxdt = strcmpi(dstest.allDrugs, drug);
        aat = dstest.AAMat(idxdt, :)';
        idxdcv = find(strcmpi(cvcclethaa.DRUG, drug));
        thaas0 = cvcclethaa.THLO(idxdcv);
        thaas1 = cvcclethaa.THHI(idxdcv);
%         for j=1:length(thaas0)
            neglabel = aat <= thaas0; 
            poslabel = aat >= thaas1;
            tp = sum(classtest == 1 & poslabel == 1);
            fp = sum(classtest == 1 & neglabel == 1);
            tn = sum(classtest == 0 & neglabel == 1);
            fn = sum(classtest == 0 & poslabel == 1);
            fprintf(fout, ',%d,%d,%d,%d', tp, fp, tn, fn);
%         end
            accdiff(end+1) = (tp+tn)/(tp+fp+tn+fn);
    catch err
        %fprintf('%s: %s', drug, err.message);
%         for j=1:length(aacutoffs)
            fprintf(fout, ',NaN,NaN,NaN,NaN');
            accdiff(end+1) = nan;
%         end
    end
    fprintf(fout, '\n');
%     fprintf('\n');
    
    %fprintf('%s, NOFS = %.2f, DIFF = %.2f\n', drug, accnofs(i), accdiff(i));
end

idxnan = isnan(accnofs) | isnan(accfs) | isnan(accdiff);

accnofs = accnofs(~idxnan);
accfs = accfs(~idxnan);
accdiff = accdiff(~idxnan);
fprintf('Num test = %d/%d\n', sum(~idxnan), length(idxnan));
fprintf('DIFF > NOFS = %d, DIFF | FS > NOFS = %d, DIFF < NOFS = %d\n', sum(accnofs<accdiff), sum(accfs > accnofs | accdiff > accnofs), sum(accnofs > accdiff));


[h, p] = ttest(accdiff, accnofs);

fprintf('mean(NOFS) = %.2f, mean(DIFF) = %.2f, p-val = %.2f\n', mean(accnofs), mean(accdiff), p);
fclose(fout);