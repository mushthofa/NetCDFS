% load dataload
% load kfgdscnet2.mat
% seldrugs = readtable('data/seldrugs.csv', 'Delimiter', ',');
% gdsc_ic50th = readtable('data/GDSC_TH.csv');
% maxn = 3;
% thaa = [0.01 0.4];
% kfold = 5;


%confmattbl = struct(); %table('', 0, 0, 0, 0, 0, 0, 'VariableNames', ...

% %{'DRUG_NAME', 'TP', 'TN', 'FP', 'FN', 'SENS', 'SPEC'});
% DRUG_NAME = {};
% TP = [];
% TN= [];
% FP= [];
% FN = [];
% SPEC = [];
% SENS = [];
% 
% for i=1:length(gdscnet.allDrugs)
%     drug = gdscnet.allDrugs{i};
%     if(sum(strcmpi(seldrugs.DRUG_NAME, drug)>0))
%         drug
%         tbl = topGenes(gdsc_dgnet{i}, gdscnet, net, 0.05, 1e-5);
%         writetable(tbl, sprintf('newlists/kfgdscnet2_%s.csv', gdscnet.allDrugs{i}));
%         idxth = strcmpi(gdsc_ic50th.Drug_Name, gdscnet.allDrugs{i}); 
%         if(sum(idxth)==0)
%             continue;
%         end
%         thic50 = gdsc_ic50th.Threshold(idxth);
%         thic50 = repmat(thic50, 1, 2);
%         [pmat, rowlabels, collabels, confmat, predImp, treemdl] = classifyCV(gdscnet, drug, tbl, maxn, thaa, thic50, kfold);
%         view(treemdl, 'mode', 'graph');
%         print(sprintf('newlists/%s_tree.png', drug), '-dpng');
%         DRUG_NAME{end+1} = drug;
%         TP(end+1) = confmat.tp;
%         TN(end+1) = confmat.tn;
%         FP(end+1) = confmat.fp;
%         FN(end+1) = confmat.fn;
%         SENS(end+1) = confmat.tp/(confmat.tp + confmat.tn);
%         SPEC(end+1) = confmat.tn/(confmat.tn + confmat.fp);
%     end
% end
% confmattbl = table(DRUG_NAME', TP', TN', FP', FN', SENS', SPEC', 'VariableNames', {'DRUG_NAME'; 'TP';'TN'; 'FP';'FN';'SENS';'SPEC'});
% writetable(confmattbl, 'newlists/kfgdscnet2_confmat.csv');


load allgdscnet1.mat
maxn = 3;
thaa = [0.05 0.4];
thic50 = [1 -2];


DRUG_NAME = {};
TP = [];
TN= [];
FP= [];
FN = [];
SPEC = [];
PREC = [];
REC = [];
for i=1:length(gdscnet.allDrugs)
    drug = gdscnet.allDrugs{i};
    if(sum(strcmpi(seldrugs.DRUG_NAME, drug)>0))
        drug
        tbl = topGenes(gdsc_dgnet{i}, gdscnet, net, 0.05, 1e-5);
        try
            [pmat, rowlabels, collabels, confmat, predImp, treemdl, predlabel1, selmat, sellabels] = classifyCV(gdscnet, drug, tbl, maxn, thaa, thic50, kfold);
        catch err
            continue;
        end
        writetable(tbl, sprintf('newlists2/ranking_%s.csv', gdscnet.allDrugs{i}));
        hg = HeatMap(selmat, 'RowLabels', sellabels);
        hg.addTitle(sprintf('Predictor matrix for %s', drug));
        ax = hg.plot;
        colorbar('peer', ax);
        fig = gcf;
        fig.Position(3:4) = 1.5*fig.Position(3:4);
        saveas(gcf, sprintf('newlists2/%s_predmat.png',drug));
        close;
        
        aa = gdscnet.AAMat(i, :);
        predlabel = predlabel1;
        predlabel(isnan(aa)) = [];
        aa(isnan(aa)) = [];
        aap = aa(strcmp(predlabel, 'p'));
        aan = aa(strcmp(predlabel, 'n'));
        barx = 0:0.05:1;
        freqp = histc(aap, barx);
        freqn = histc(aan, barx);
        bar(barx, [freqp' freqn']);
        set(gca, 'XLim', [0 1])
        legend('Positive', 'Negative');
        title(sprintf('Prediction for %s', drug));
        xlabel('Activity Area');
        ylabel('Frequency');
        print(sprintf('newlists2/%s_aapred.png', drug), '-dpng');
        close;
        
        predlabel = predlabel1;
        ic50 = gdscnet.IC50Mat(i, :);
        predlabel(isnan(ic50)) = [];
        ic50(isnan(ic50)) = [];
        ic50p = ic50(strcmp(predlabel, 'p'));
        ic50n = ic50(strcmp(predlabel, 'n'));
        mini = min(ic50);
        maxi = max(ic50);
        barx = mini:0.5:maxi;
        freqp = histc(ic50p, barx);
        freqn = histc(ic50n, barx);
        bar(barx, [freqp' freqn']);
        set(gca, 'XLim', [mini maxi])
        legend('Positive', 'Negative');
        title(sprintf('Prediction for %s', drug));
        xlabel('IC50');
        ylabel('Frequency');
        print(sprintf('newlists2/%s_ic50pred.png', drug), '-dpng');
        close;
        
        h = treemdl.Impl.viewGraph([], strcmp(treemdl.NodeClass, 'p'),treemdl.PredictorNames,prunelevel,'');
        set(h, 'Position', [0 0 800 600]);
        saveas(h, sprintf('newlists2/%s_tree.png', drug));
        close(h);
        %print(sprintf('newlists/%s_tree.png', drug), '-dpng');
        DRUG_NAME{end+1} = drug;
        TP(end+1) = confmat.tp;
        TN(end+1) = confmat.tn;
        FP(end+1) = confmat.fp;
        FN(end+1) = confmat.fn;
        SPEC(end+1) = confmat.tn/(confmat.tn + confmat.fp);
        PREC(end+1) = confmat.tp/(confmat.tp + confmat.fp);
        REC(end+1) = confmat.tp/(confmat.tp + confmat.fn);
    end
end
confmattbl = table(DRUG_NAME', TP', TN', FP', FN', SPEC', PREC', REC', 'VariableNames', {'DRUG_NAME'; 'TP';'TN'; 'FP';'FN';'SPEC'; 'PREC'; 'REC'});
writetable(confmattbl, 'newlists2/allgdscnet1_confmat.csv');





