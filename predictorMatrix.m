function [pmat, mutfeats, cnvfeats, gexfeats] = predictorMatrix(ds, feats)

pmat = [];
cnvfeats = cell(0);
mutfeats = cell(0);
gexfeats = cell(0);
for i=1:length(feats)
    
    ff = strsplit(feats{i}, '-');
    genename = strjoin(ff(1:end-1), '-');
    ptype = ff{end};
    switch(ptype)
        case 'MUT'
            idxg = strcmpi(ds.mutGenes, genename);
            if(sum(idxg)~=1)
                error('Feature %s not found\n', feats{i});
                pmat = [pmat; zeros(1, length(ds.cellNames))];
            else
                pmat = [pmat; ds.mutMat(idxg, :)];
                mutfeats = [mutfeats; genename];
            end
        case 'CNV'
            idxg = strcmpi(ds.cnvGenes, genename);
            if(sum(idxg)~=1)
                error('Feature %s not found\n', feats{i});
                pmat = [pmat; zeros(1, length(ds.cellNames))];
            else
                pmat = [pmat; ds.cnvMat(idxg, :)];
                cnvfeats = [cnvfeats; genename];
            end
        case 'GEX'
            idxg = strcmpi(ds.dgexGenes, genename);
            if(sum(idxg)~=1)
                error('Feature %s not found\n', feats{i});
                pmat = [pmat; zeros(1, length(ds.cellNames))];
            else
                pmat = [pmat; ds.dgexMat(idxg, :)];
                gexfeats = [gexfeats; genename];
            end
        otherwise
            pmat = NaN;
            error('Unknown predictor type %s', ptype);
    end
    
end

end