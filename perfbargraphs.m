% Bargraphs for performances on GDSC CV

% Overall performance/accuracy
nofsacc = (cvtab.TPNOFS + cvtab.TNNOFS)./(cvtab.TPNOFS + cvtab.TNNOFS + cvtab.FPNOFS + cvtab.FNNOFS);
fsacc = (cvtab.TPFS + cvtab.TNFS)./(cvtab.TPFS + cvtab.TNFS + cvtab.FPFS + cvtab.FNFS);
diffacc = (cvtab.TPDIFF + cvtab.TNDIFF)./(cvtab.TPDIFF + cvtab.TNDIFF + cvtab.FPDIFF + cvtab.FNDIFF);
bar([nofsacc, fsacc, diffacc]);
legend('NOFS', 'FS', 'FS+Diff')
set(gca, 'XTickLabel', cvtab.DRUG)
set(gca, 'XTick', 1:height(cvtab))
xtickangle(45)
title('Accuracy of models on GDSC cross validations');
ylabel('Accuracy');
bpmsg1 = sprintf('NOFS (mean=%.3f)', mean(nofsacc));
bpmsg2 = sprintf('FS (mean=%.3f)', mean(fsacc));
bpmsg3 = sprintf('FS + Diff (mean=%.3f)', mean(diffacc));
boxplot([nofsacc, fsacc, diffacc], {bpmsg1, bpmsg2, bpmsg3});
[~, p] = ttest(nofsacc, diffacc);
title(sprintf('Accuracy distribution, T-test p-val(NOFS, FS+DIFF) = %0.5e', p));
ylabel('Accuracy');

% Positive recalls
nofsposrec = (cvtab.TPNOFS)./(cvtab.TPNOFS + cvtab.FNNOFS);
fsposrec = (cvtab.TPFS)./(cvtab.TPFS + cvtab.FNFS);
diffposrec = (cvtab.TPDIFF)./(cvtab.TPDIFF + cvtab.FNDIFF);
bar([nofsposrec, fsposrec, diffposrec]);
legend('NOFS', 'FS', 'FS+Diff')
set(gca, 'XTickLabel', cvtab.DRUG)
set(gca, 'XTick', 1:height(cvtab))
xtickangle(45);
title('Recall on positive response on GDSC cross validations');
ylabel('Recall');
bpmsg1 = sprintf('NOFS (mean=%.2f)', mean(nofsposrec));
bpmsg2 = sprintf('FS (mean=%.2f)', mean(fsposrec));
bpmsg3 = sprintf('FS + Diff (mean=%.2f)', mean(diffposrec));
boxplot([nofsposrec, fsposrec, diffposrec], {bpmsg1, bpmsg2, bpmsg3});
[h, p] = ttest(nofsposrec, diffposrec);
title(sprintf('Recall on pos. response, T-test p-value between NOFS and FS+DIFF = %0.5e', p));
ylabel('Recall');

% Negative recalls
nofsnegrec = (cvtab.TNNOFS)./(cvtab.TNNOFS + cvtab.FPNOFS);
fsnegrec = (cvtab.TNFS)./(cvtab.TNFS + cvtab.FPFS);
diffnegrec = (cvtab.TNDIFF)./(cvtab.TNDIFF + cvtab.FPDIFF);
bar([nofsnegrec, fsnegrec, diffnegrec]);
legend('NOFS', 'FS', 'FS+Diff')
set(gca, 'XTickLabel', cvtab.DRUG)
set(gca, 'XTick', 1:height(cvtab))
xtickangle(45);
title('Accuracy of models on GDSC Cross Validations');
ylabel('Accuracy');
bpmsg1 = sprintf('NOFS (mean=%.2f)', mean(nofsnegrec));
bpmsg2 = sprintf('FS (mean=%.2f)', mean(fsnegrec));
bpmsg3 = sprintf('FS + Diff (mean=%.2f)', mean(diffnegrec));
boxplot([nofsnegrec, fsnegrec, diffnegrec], {bpmsg1, bpmsg2, bpmsg3});
[h, p] = ttest(nofsnegrec, diffnegrec);
title(sprintf('Recall on neg. response, T-test p-value between NOFS and FS+DIFF = %0.5e', p));
ylabel('Recall');
