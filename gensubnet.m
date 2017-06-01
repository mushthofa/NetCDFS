rng('default');
minleafsize = 5;
eps = 1e-5;
maxiter = 30;
%maxn = 10;
nrep = 1;
%ntrees = 30;
kfold = 5;
ds = gdscnet;
allfeats = readtable('manth/allfeats.csv');
kflnofs = zeros(height(seldrugs), 1);
kflfs = kflnofs;
kfldiff = kflnofs;
fout = fopen(sprintf('subnets/diffkfold_rng_def_%.2f_%d_%d.csv', alpha, maxn, ntrees), 'w');
fprintf(fout, 'No,DRUG,FRACN,FRACP,TPNOFS,FPNOFS,TNNOFS,FNNOFS,TPFS,FPFS,TNFS,FNFS,TPDIFF,FPDIFF,TNDIFF,FNDIFF\n');
for i=1:height(seldrugs)
    drug = seldrugs.DRUG_NAME{i};
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
    
    
    % Train classification tree without feature selection
    
    featmat = [];
    featmat = [featmat; ds.mutMat(:, [find(idxneg) find(idxpos)])];
    featmat = [featmat; ds.cnvMat(:, [find(idxneg) find(idxpos)])];
    featmat = [featmat; ds.dgexMat(:, [find(idxneg) find(idxpos)])];
    alltissues = unique(ds.cellTissues);
    tissmat = zeros(length(alltissues), length(ds.cellNames));
    for j=1:length(alltissues)
        tt = alltissues{j};
        tissmat(j, :) = strcmpi(ds.cellNames, tt)';
    end
    featmat = [featmat; tissmat(:, [find(idxneg) find(idxpos)])];
    featmat = featmat';
    classtarget = [zeros(sum(idxneg), 1); ones(sum(idxpos), 1)];
    pnames = strcat(ds.mutGenes, '-MUT');
    pnames = [pnames; strcat(ds.cnvGenes, '-CNV')];
    pnames = [pnames; strcat(ds.dgexGenes, '-GEX')];
    pnames = [pnames; strcat('tissue_is_', alltissues)];
    
    tpnofs = zeros(nrep, 1);
    fpnofs = zeros(nrep, 1);
    tnnofs = zeros(nrep, 1);
    fnnofs = zeros(nrep, 1);
    
    for j=1:nrep
        tt = templateTree('MinLeafSize', minleafsize);
        treemdl = fitcensemble(featmat, classtarget, 'Learners', tt, 'Method', 'GentleBoost', 'NumLearningCycles', ntrees, 'CrossVal', 'on', 'KFold', kfold, 'PredictorNames', pnames);
        label = kfoldPredict(treemdl);
        tpnofs(j) = sum(classtarget == 1 & label == 1);
        fpnofs(j) = sum(classtarget == 0 & label == 1);
        tnnofs(j) = sum(classtarget == 0 & label == 0);
        fnnofs(j) = sum(classtarget == 1 & label == 0);
    end
    
    
    
    % Train classification tree WITH feature selection
    
    idxmut = ismember(ds.mutGenes, allgenes);
    idxcnv = ismember(ds.cnvGenes, allgenes);
    idxgex = ismember(ds.dgexGenes, allgenes);
    
    featmat = [];
    featmat = [featmat; ds.mutMat(idxmut, [find(idxneg) find(idxpos)])];
    featmat = [featmat; ds.cnvMat(idxcnv, [find(idxneg) find(idxpos)])];
    featmat = [featmat; ds.dgexMat(idxgex, [find(idxneg) find(idxpos)])];
    
    featmat = [featmat; tissmat(:, [find(idxneg) find(idxpos)])];
    featmat = featmat';
    classtarget = [zeros(sum(idxneg), 1); ones(sum(idxpos), 1)];
    pnames = strcat(ds.mutGenes(idxmut), '-MUT');
    pnames = [pnames; strcat(ds.cnvGenes(idxcnv), '-CNV')];
    pnames = [pnames; strcat(ds.dgexGenes(idxgex), '-GEX')];
    pnames = [pnames; strcat('tissue_is_', alltissues)];
    
    tpfs = zeros(nrep, 1);
    fpfs = zeros(nrep, 1);
    tnfs = zeros(nrep, 1);
    fnfs = zeros(nrep, 1);
    
    for j=1:nrep
        tt = templateTree('MinLeafSize', minleafsize);
        treemdl = fitcensemble(featmat, classtarget, 'Learners', tt, 'Method', 'GentleBoost', 'NumLearningCycles', ntrees, 'CrossVal', 'on', 'KFold', kfold, 'PredictorNames', pnames);
        label = kfoldPredict(treemdl);
        tpfs(j) = sum(classtarget == 1 & label == 1);
        fpfs(j) = sum(classtarget == 0 & label == 1);
        tnfs(j) = sum(classtarget == 0 & label == 0);
        fnfs(j) = sum(classtarget == 1 & label == 0);
    end
    
    % Classification tree but diffuse mutation + CNV first
    

    % Diffuse mutations
    mutgenes = union(mutpgenes, mutngenes);
    idxmut = cellfun(@(x) find(ismember(ds.mutGenes, x)), mutgenes);
    mutmat = ds.mutMat(idxmut, [find(idxneg) find(idxpos)]);
    
    
    idxmutnet = cell2mat(values(net.name2node, mutgenes));
    adjmat = net.mat(idxmutnet, idxmutnet);
    adjmat = single((adjmat + adjmat')>0);
    
    nSample = sum(idxpos)+sum(idxneg);
    diffmutMat = netSmooth(mutmat', adjmat, alpha, eps, maxiter)';
    diffmutMat = double(diffmutMat >= (1-alpha));
    
    
    % Diffuse Amplifications
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
    
    
    featmat = [featmat; tissmat(:, [find(idxneg) find(idxpos)])];
    featmat = featmat';
    
%     Zero/randomise
%     featmat = zeros(size(featmat));
    
    
    pnames = strcat(mutgenes, '-MUT');
    pnames = [pnames; strcat(ds.cnvGenes(idxcnv), '-CNV')];
    pnames = [pnames; strcat(ds.dgexGenes(idxgex), '-GEX')];
    pnames = [pnames; strcat('tissue_is_', alltissues)];
    
    tpdiff = zeros(nrep, 1);
    fpdiff = zeros(nrep, 1);
    tndiff = zeros(nrep, 1);
    fndiff = zeros(nrep, 1);
    
    for j=1:nrep
        tt = templateTree('MinLeafSize', minleafsize);
        treemdl = fitcensemble(featmat, classtarget, 'Learners', tt, 'Method', 'GentleBoost', 'NumLearningCycles', ntrees, 'CrossVal', 'on', 'KFold', kfold, 'PredictorNames', pnames);
        label = kfoldPredict(treemdl);
        tpdiff(j) = sum(classtarget == 1 & label == 1);
        fpdiff(j) = sum(classtarget == 0 & label == 1);
        tndiff(j) = sum(classtarget == 0 & label == 0);
        fndiff(j) = sum(classtarget == 1 & label == 0);
    end;
    
    
    fprintf(fout, '%d,%s,%.2f,%.2f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n', ...
        i, drug, sum(idxneg)/(sum(idxneg)+sum(idxpos)), sum(idxpos)/(sum(idxneg)+sum(idxpos)), mean(tpnofs), mean(fpnofs), mean(tnnofs), mean(fnnofs),  mean(tpfs), mean(fpfs), mean(tnfs), mean(fnfs), mean(tpdiff), mean(fpdiff), mean(tndiff), mean(fndiff));
    fprintf('%d,%s,%.2f,%.2f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n', ...
        i, drug, sum(idxneg)/(sum(idxneg)+sum(idxpos)), sum(idxpos)/(sum(idxneg)+sum(idxpos)), mean(tpnofs), mean(fpnofs), mean(tnnofs), mean(fnnofs),  mean(tpfs), mean(fpfs), mean(tnfs), mean(fnfs), mean(tpdiff), mean(fpdiff), mean(tndiff), mean(fndiff))
end
fclose(fout);


