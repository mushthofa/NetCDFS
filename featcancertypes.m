allfeatnames = readtable('manth/allfeats.csv', 'Delimiter', ',');
seldrugs = readtable('data/seldrugs.csv', 'Delimiter', ',');
fout = fopen('manth/featcancertypes.csv', 'w');
fprintf(fout, 'DRUG, FEATURE, SIGN, TP, TN, FP, FN, RECALL, PREC, Cancer Types\n');
for i=1:height(allfeatnames)
    drug = allfeatnames.DRUG_NAME{i};
    feat = allfeatnames.FEATURE{i};
    sign = allfeatnames.SIGN{i};
    idxd = strcmpi(seldrugs.DRUG_NAME, drug);
    aa1 = seldrugs.AATH1(idxd);
    aa2 = seldrugs.AATH2(idxd);
    
    idxd = strcmpi(gdscnet.allDrugs, drug);
    idxpos = gdscnet.AAMat(idxd, :) >= aa2;
    idxneg = gdscnet.AAMat(idxd, :) <= aa1;
    
    idxsplit = find(feat=='-');
    idxsplit = idxsplit(end);
    gene = feat(1:idxsplit-1);
    type = feat(idxsplit+1:end);
    
    if(strcmp(type, 'MUT'))
        idxfeat = strcmpi(gdscnet.mutGenes, gene);
        datafeat = gdscnet.mutMat;
    elseif(strcmp(type, 'CNV') || strcmp(type, 'DEL') || strcmp(type, 'AMP'))
        idxfeat = strcmpi(gdscnet.cnvGenes, gene);
        datafeat = gdscnet.cnvMat~=0;
    elseif(strcmp(type, 'GEX'))
        idxfeat = strcmpi(gdscnet.dgexGenes, gene);
        datafeat = gdscnet.dgexMat~=0;
    end
    
    if(sign == '+')
        posset = gdscnet.cellOriTissues(idxpos & datafeat(idxfeat, :)~=0);
%         negset = gdscnet.cellOriTissues(idxneg & datafeat(idxfeat, :)==0);
        tp = sum(idxpos & datafeat(idxfeat, :)~=0);
        tn = sum(idxneg & datafeat(idxfeat, :)==0);
        fp = sum(idxpos & datafeat(idxfeat, :)==0);
        fn = sum(idxneg & datafeat(idxfeat, :)~=0);
        
    elseif (sign == '-')
        posset = gdscnet.cellOriTissues(idxneg & datafeat(idxfeat, :)~=0);
%         negset = gdscnet.cellOriTissues(idxneg & datafeat(idxfeat, :)~=0);
        tp = sum(idxneg & datafeat(idxfeat, :)~=0);
        tn = sum(idxpos & datafeat(idxfeat, :)==0);
        fp = sum(idxneg & datafeat(idxfeat, :)==0);
        fn = sum(idxpos & datafeat(idxfeat, :)~=0);
    end
    rec = tp/(tp+fn);
    prec = tp/(tp+fp);
    ttpos = tabulate(posset);
%     ttneg = tabulate(negset);
    ttpos = sortrows(ttpos, -2);
%     ttneg = sortrows(ttneg, -2);
    fprintf(fout, '%s, %s, %s, %d, %d, %d, %d, %.2f, %.2f,  ', drug, feat, sign, tp, tn, fp, fn, rec, prec);
    for j=1:size(ttpos, 1)
        fprintf(fout, '%s(%d) ', ttpos{j, 1}, ttpos{j, 2});
    end
%     fprintf(fout, ',');
%     for j=1:size(ttneg, 1)
%         fprintf(fout, '%s(%d) ', ttneg{j, 1}, ttneg{j, 2});
%     end
    fprintf(fout, '\n');
end
fclose(fout);