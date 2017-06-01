function [pmat, rowlabels, collabels, predImp, treemdl, predlabel, predprob, selmat, sellabels, rf, rfpred, rfprob] = classifyCV(ds, drug, tbl, maxn, thaa, kfold)

mutpgenes = tbl.MUTP;
mutngenes = tbl.MUTN;
cnvpgenes = tbl.CNVP;
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

idxdrug = strcmpi(ds.allDrugs, drug);
[~, idxsort] = sort(ds.AAMat(idxdrug, :));
idxrest = find(ds.AAMat(idxdrug, idxsort)<=thaa(1)); %  & ds.IC50Mat(idxdrug, idxsort)>thic50(2)); %;  | ds.IC50Mat(idxdrug, idxsort)>=thic50(1));
idxresp = find(ds.AAMat(idxdrug, idxsort)>=thaa(2)); % | ds.IC50Mat(idxdrug, idxsort)<=thic50(2));

if(isempty(idxrest) || isempty(idxresp))
    error('Error... empty samples in drug %s ', drug);
end
idxmutp = cell2mat(cellfun(@(x) find(strcmp(x, ds.mutGenes)), mutpgenes, 'UniformOutput', 0));
idxcnvp = cell2mat(cellfun(@(x) find(strcmp(x, ds.cnvGenes)), cnvpgenes, 'UniformOutput', 0));
idxgexp = cell2mat(cellfun(@(x) find(strcmp(x, ds.dgexGenes)), gexpgenes, 'UniformOutput', 0));
idxmutn = cell2mat(cellfun(@(x) find(strcmp(x, ds.mutGenes)), mutngenes, 'UniformOutput', 0));
idxcnvn = cell2mat(cellfun(@(x) find(strcmp(x, ds.cnvGenes)), cnvngenes, 'UniformOutput', 0));
idxgexn = cell2mat(cellfun(@(x) find(strcmp(x, ds.dgexGenes)), gexngenes, 'UniformOutput', 0));

pmat = []; %ds.cellTissues(idxsort([idxrest idxresp]));
pmat = [pmat; ds.mutMat([idxmutp; idxmutn], idxsort([idxrest idxresp]))];
pmat = [pmat; ds.cnvMat([idxcnvp; idxcnvn], idxsort([idxrest idxresp]))];
pmat = [pmat; ds.dgexMat([idxgexp; idxgexn], idxsort([idxrest idxresp]))];


% Add tissue labels
alltissues = unique(ds.cellTissues);
tissmat = zeros(length(alltissues), length(ds.cellNames));
for i=1:length(alltissues)
    tt = alltissues{i};
    tissmat(i, :) = strcmpi(ds.cellNames, tt)';
end
pmat = [pmat; tissmat(:, idxsort([idxrest idxresp]))];
sampleIn = pmat';


classOut = repmat({'n'}, length(idxrest), 1);
classOut = [classOut; repmat({'p'}, length(idxresp), 1)];




% tp = 0;
% tn = tp;
% fp = tp;
% fn = tp;
% % predImp = zeros(size(sampleIn, 2), kfold);
% % for i=1:kfold
% %     selind = crossvalind('KFold', nData, kfold);
% %     testIdx = selind==i;
% %     trainIdx = ~testIdx;
% %     tr = fitensemble(sampleIn(trainIdx, :), classOut(trainIdx), 'Bag', ntrees, 'Tree', 'Type', 'C');
% %     predOut = tr.predict(sampleIn(testIdx, :));
% %     testOut = classOut(testIdx);
% %     tp = tp + sum(strcmp(testOut, 'p') & strcmp(predOut, 'p'));
% %     tn = tn + sum(strcmp(testOut, 'n') & strcmp(predOut, 'n'));
% %     fp = fp + sum(strcmp(testOut, 'n') & strcmp(predOut, 'p'));
% %     fn = fn + sum(strcmp(testOut, 'p') & strcmp(predOut, 'n'));
% %     predImp(:, i) = predictorImportance(tr);
% % end
% % predImp = sum(predImp, 2);
% confmat.tp = tp;
% confmat.tn = tn;
% confmat.fp = fp;
% confmat.fn = fn;



% 
ntrees = 500;
rowlabels = strcat([mutpgenes; mutngenes], '-MUT');
rowlabels = [rowlabels; strcat([cnvpgenes; cnvngenes], '-CNV')];
rowlabels = [rowlabels; strcat([gexpgenes; gexngenes], '-GEX')];
rowlabels = [rowlabels; strcat('tissue_is_', alltissues)];

numsmallclass = min(sum(idxrest), sum(idxresp));
minleafsize = 5;  %min(floor(0.25*numsmallclass), 10);
tt = templateTree('MinLeafSize', minleafsize);
rf = fitcensemble(sampleIn, classOut, 'Method', 'Bag', 'NumLearningCycles', ntrees, 'Learners', tt); %TreeBagger(ntrees, sampleIn, classOut, 'Method', 'Classification', 'MinLeafSize', 5, 'CategoricalPredictors', 'all', 'PredictorNames', rowlabels);
predImp = predictorImportance(rf);

zpi = zscore(predImp);
pv = 1 - normcdf(abs(zpi), 0, 1);

ntake = 5; %max(sum(pv<=0.01), 5);

[~, idxsort] = sort(predImp, 'descend');
treemdl = fitctree(sampleIn(:, idxsort(1:ntake)), classOut, 'PredictorNames', rowlabels(idxsort(1:ntake)), 'MaxNumCategories', 2, 'MinLeafSize', 5, 'MaxNumSplits', 5, 'CategoricalPredictors', 'all');
sellabels = rowlabels(idxsort(1:ntake));

allmat = []; %ds.cellTissues(idxsort([idxrest idxresp]));
allmat = [allmat; ds.mutMat([idxmutp; idxmutn], :)];
allmat = [allmat; ds.cnvMat([idxcnvp; idxcnvn], :)];    
allmat = [allmat; ds.dgexMat([idxgexp; idxgexn],:)];
allmat = [allmat; tissmat];
[rfpred, rfscores] = predict(rf, allmat');
rfprob = max(rfscores')';

allmat = allmat(idxsort(1:ntake), :)';

[predlabel, cp] = predict(treemdl, allmat);
predprob = zeros(length(predlabel), 1);
predprob(strcmp(predlabel, 'n')) = cp(strcmp(predlabel, 'n'), 1);
predprob(strcmp(predlabel, 'p')) = cp(strcmp(predlabel, 'p'), 2);

aa = ds.AAMat(idxdrug, :);
allmat(isnan(aa), :) = [];
aa(isnan(aa)) = [];
[~, idxsort] = sort(aa);
allmat = allmat(idxsort, :);
allmat = allmat(:, end:-1:1);
sellabels = sellabels(end:-1:1);
selmat = allmat';

collabels = classOut;


