%load cur1
manpval = readtable('manth/manth_pval.csv');
sharedcorr = readtable('sharedcellcorr.csv') ;
idxcv = ~isnan(manpval.CTRP_AUC);
nf = 1; %sum(idxcv);
thpval = 0.05;
bluedrugs = manpval.DRUG_NAME(manpval.CTRP_AUC <= thpval/nf); % Bonferroni correction;
reddrugs = manpval.DRUG_NAME(manpval.CTRP_AUC > thpval/nf);
idxblue = ismember(sharedcorr.DRUG, bluedrugs);
idxred = ismember(sharedcorr.DRUG, reddrugs);
hold on

plot(sharedcorr.N2(idxblue), sharedcorr.CORR2(idxblue), 'ob');
plot(sharedcorr.N2(idxred), sharedcorr.CORR2(idxred), 'xr');
text(sharedcorr.N2(idxblue) + 5, sharedcorr.CORR2(idxblue), sharedcorr.DRUG(idxblue))
text(sharedcorr.N2(idxred) + 5, sharedcorr.CORR2(idxred), sharedcorr.DRUG(idxred));
xlabel('Number of shared cell lines between GDSC - CTRP');
ylabel('Drug profile correlation on shared cell lines');
title('Distribution of shared cell lines and their profile correlation between GDSC - CTRP');
legend('Significant CV P-value', 'Non-significant CV P-value');