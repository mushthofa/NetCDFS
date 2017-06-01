for i=1:height(seldrugs)
    drug = seldrugs.DRUG_NAME{i};
    idxd1 = strcmpi(gdscnet.allDrugs, drug);
    cell_gdsc = gdscnet.cellNames(~isnan(gdscnet.AAMat(idxd1, :)));
    idxd2 = strcmpi(cclenet.allDrugs, drug);
    cell_ccle = cclenet.cellNames(~isnan(cclenet.AAMat(idxd2, :)));
    idxd3 = strcmpi(ctrpnet.allDrugs, drug);
    cell_ctrp = ctrpnet.cellNames(~isnan(ctrpnet.AAMat(idxd3, :)));
    sharedcells  = intersect(cell_gdsc, cell_ctrp);
    
    idxc = cellfun(@(x) find(strcmpi(gdscnet.cellNames, x)), sharedcells);
    AA_gdsc = gdscnet.AAMat(idxd1, idxc)';
    
    idxc = cellfun(@(x) find(strcmpi(cclenet.cellNames, x)), sharedcells);
    AA_ccle = cclenet.AAMat(idxd2, idxc)';
    
    idxc = cellfun(@(x) find(strcmpi(ctrpnet.cellNames, x)), sharedcells);
    AA_ctrp = ctrpnet.AAMat(idxd3, idxc)';
    
    idx_gdsc = ~isnan(AA_gdsc);
    idx_ccle = ~isnan(AA_ccle);
    idx_ctrp = ~isnan(AA_ctrp);
    
    n1 = sum(idx_gdsc & idx_ccle);
    if(~isempty(n1) & n1 > 0)
        cc1 = corr(AA_gdsc(idx_gdsc & idx_ccle), AA_ccle(idx_gdsc & idx_ccle));
    else
        cc1 = nan;
        n1 = 0;
    end
    
    n2 = sum(idx_gdsc & idx_ctrp);
    if(~isempty(n2) & n2>0)
        cc2 = corr(AA_gdsc(idx_gdsc & idx_ctrp), AA_ctrp(idx_gdsc & idx_ctrp));
    else
        cc2 = nan;
        n2 = 0;
    end
    
    fprintf('%s, %d, %.2f, %d, %.2f\n', drug, n1, cc1, n2, cc2);
    
end