corfac = 1; %sum(~isnan(pvals.CCLE_AUC));
idxsel = ~isnan(pvals.CCLE_AUC);
[sortedpvals, idxsort] = sort(-log10(pvals.CCLE_AUC(idxsel)));
sortedpvals(sortedpvals>30) = 30;
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
thpval = 0.05;
disp('CCLE');
numred = sum(sortedpvals < -log10(thpval/corfac));
numblue = sum(sortedpvals >  -log10(thpval/corfac));
fprintf('Validated = %d/%d\n', numblue, numred+numblue);
h1 = bar(1:numred, sortedpvals(sortedpvals < -log10(thpval/corfac)), 'parent', aHand, 'Facecolor', 'r');
h2 = bar(numred+1:numred+numblue, sortedpvals(sortedpvals >= -log10(thpval/corfac)), 'parent', aHand, 'Facecolor', 'b');
labels = pvals.DRUG_NAME(idxsel);
set(gca, 'XTickLabel', labels(idxsort))
set(gca, 'XTick', 1:length(sortedpvals))
set(gca, 'XTickLabelRotation', 45);
ylabel('-log_{10}(p-value)');
title('Rank sum test for validation on CCLE');
legend('P-value > 0.05', 'P-value <=0.05');

corfac = 1; %sum(~isnan(pvals.CCLE_AUC));
idxsel = ~isnan(pvals.CTRP_AUC);
[sortedpvals, idxsort] = sort(-log10(pvals.CTRP_AUC(idxsel)));
sortedpvals(sortedpvals>30) = 30;
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
thpval = 0.05;
numred = sum(sortedpvals < -log10(thpval/corfac));
numblue = sum(sortedpvals >  -log10(thpval/corfac));
disp('CTRP');
fprintf('Validated = %d/%d\n', numblue, numred+numblue);

h1 = bar(1:numred, sortedpvals(sortedpvals < -log10(thpval/corfac)), 'parent', aHand, 'Facecolor', 'r');
h2 = bar(numred+1:numred+numblue, sortedpvals(sortedpvals >= -log10(thpval/corfac)), 'parent', aHand, 'Facecolor', 'b');
labels = pvals.DRUG_NAME(idxsel);
set(gca, 'XTickLabel', labels(idxsort))
set(gca, 'XTick', 1:length(sortedpvals))
set(gca, 'XTickLabelRotation', 45);
ylabel('-log_{10}(p-value)');
title('Rank sum test for validation on CTRP');
legend('P-value > 0.05', 'P-value <=0.05');

corfac = 1; %sum(~isnan(pvals.CCLE_AUC));
idxsel = ~isnan(pvals.GDSC_AUC);
[sortedpvals, idxsort] = sort(-log10(pvals.GDSC_AUC(idxsel)));
sortedpvals(sortedpvals>30) = 30;
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
thpval = 0.05;
numred = sum(sortedpvals < -log10(thpval/corfac));
numblue = sum(sortedpvals >  -log10(thpval/corfac));
disp('GDSC');
fprintf('Validated = %d/%d\n', numblue, numred+numblue);
h1 = bar(1:numred, sortedpvals(sortedpvals < -log10(thpval/corfac)), 'parent', aHand, 'Facecolor', 'r');
h2 = bar(numred+1:numred+numblue, sortedpvals(sortedpvals >= -log10(thpval/corfac)), 'parent', aHand, 'Facecolor', 'b');
labels = pvals.DRUG_NAME(idxsel);
set(gca, 'XTickLabel', labels(idxsort))
set(gca, 'XTick', 1:length(sortedpvals))
set(gca, 'XTickLabelRotation', 45);
ylabel('-log_{10}(p-value)');
title('Rank sum test for validation on GDSC');
legend('P-value > 0.05', 'P-value <=0.05');