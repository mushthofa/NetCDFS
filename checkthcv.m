drugshared = intersect(ctrpnet.allDrugs, seldrugs.DRUG_NAME);
fout = fopen('cvctrpthaa.csv', 'w');
fprintf(fout, 'DRUG,THLO,THHI,CORR\n');
for i=1:length(drugshared)
    maxacc = 0;
    maxth = nan;
    for aacutoffs = 0.01:0.01:0.99
        drug = drugshared{i};
        idxd = strcmpi(ctrpnet.allDrugs, drug);
        aactrp = ctrpnet.AAMat(idxd, :);
        curth = quantile(aactrp, aacutoffs);
        poscells_ctrp = ctrpnet.cellNames(aactrp>=curth);
        negcells_ctrp = ctrpnet.cellNames(aactrp<curth);
        idxseld = strcmpi(seldrugs.DRUG_NAME, drug);
        idxd = strcmpi(gdscnet.allDrugs, drug);
        aagdsc = gdscnet.AAMat(idxd,:);
        poscells_gdsc = gdscnet.cellNames(aagdsc>=seldrugs.AATH2(idxseld));
        negcells_gdsc = gdscnet.cellNames(aagdsc<=seldrugs.AATH1(idxseld));
        tp = length(intersect(poscells_gdsc, poscells_ctrp));
        tn = length(intersect(negcells_ctrp, negcells_gdsc));
        cp =  tp + tn;
        all = length(intersect(union(poscells_ctrp, negcells_ctrp), union(poscells_gdsc, negcells_gdsc)));
        acc = cp/all;
        if(acc > maxacc)
            maxacc = acc;
            maxth = curth;
        elseif (acc==maxacc)
            maxth(end+1) = curth;
        end
    end
    fprintf(fout, '%s,%.2f,%.2f,%.2f\n', drug, median(maxth), median(maxth), maxacc);
end
fclose(fout);