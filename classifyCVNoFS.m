function [pmat, rowlabels, collabels, predImp, treemdl, predlabel, predprob, selmat, sellabels, rf, rfpred, rfprob] = classifyCVNoFS(ds, drug, thaa)

idxdrug = strcmpi(ds.allDrugs, drug);
[~, idxsort] = sort(ds.AAMat(idxdrug, :));
idxrest = find(ds.AAMat(idxdrug, idxsort)<=thaa(1)); %  & ds.IC50Mat(idxdrug, idxsort)>thic50(2)); %;  | ds.IC50Mat(idxdrug, idxsort)>=thic50(1));
idxresp = find(ds.AAMat(idxdrug, idxsort)>=thaa(2)); % | ds.IC50Mat(idxdrug, idxsort)<=thic50(2));

if(isempty(idxrest) || isempty(idxresp))
    error('Error... empty samples in drug %s ', drug);
end


pmat = []; %ds.cellTissues(idxsort([idxrest idxresp]));
pmat = [pmat; ds.mutMat(:, idxsort([idxrest idxresp]))];
pmat = [pmat; ds.cnvMat(:, idxsort([idxrest idxresp]))];
pmat = [pmat; ds.dgexMat(:, idxsort([idxrest idxresp]))];


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
rowlabels = strcat(ds.mutGenes, '-MUT');
rowlabels = [rowlabels; strcat(ds.cnvGenes, '-CNV')];
rowlabels = [rowlabels; strcat(ds.dgexGenes, '-GEX')];
rowlabels = [rowlabels; strcat('tissue_is_', alltissues)];

% numsmallclass = min(sum(idxrest), sum(idxresp));
disp('Training random forest...');
minleafsize = 5;  %min(floor(0.25*numsmallclass), 10);
tt = templateTree('MinLeafSize', minleafsize);
rf = fitcensemble(sampleIn, classOut, 'Method', 'Bag', 'NumLearningCycles', ntrees, 'Learners', tt); %TreeBagger(ntrees, sampleIn, classOut, 'Method', 'Classification', 'MinLeafSize', 5, 'CategoricalPredictors', 'all', 'PredictorNames', rowlabels);
predImp = predictorImportance(rf);

zpi = zscore(predImp);
pv = 1 - normcdf(abs(zpi), 0, 1);

ntake = 5; %max(sum(pv<=0.01), 5);

[~, idxsort] = sort(predImp, 'descend');
treemdl = fitrtree(sampleIn(:, idxsort(1:ntake)), classOut, 'PredictorNames', rowlabels(idxsort(1:ntake)), 'MaxNumCategories', 2, 'MinLeafSize', 5, 'MaxNumSplits', 5, 'CategoricalPredictors', 'all');
sellabels = rowlabels(idxsort(1:ntake));

allmat = []; %ds.cellTissues(idxsort([idxrest idxresp]));
allmat = [allmat; ds.mutMat];
allmat = [allmat; ds.cnvMat];    
allmat = [allmat; ds.dgexMat];
allmat = [allmat; tissmat];

disp('Predicting...');
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


