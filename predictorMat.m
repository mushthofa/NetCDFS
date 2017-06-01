function [pmat, rowlabels, profs] = predictorMat(ds, drug, tbl, maxn)%, thaa)%, thic50)

mutpgenes = tbl.MUTP;
mutpgenes = mutpgenes(ismember(mutpgenes, ds.mutGenes));
mutngenes = tbl.MUTN;
mutngenes = mutngenes(ismember(mutngenes, ds.mutGenes));
cnvpgenes = tbl.CNVP;
cnvpgenes = cnvpgenes(ismember(cnvpgenes, ds.cnvGenes));
cnvngenes = tbl.CNVN;
cnvngenes = cnvngenes(ismember(cnvngenes, ds.cnvGenes));
gexpgenes = tbl.GEXP;
gexpgenes = gexpgenes(ismember(gexpgenes, ds.dgexGenes));
gexngenes = tbl.GEXN;
gexngenes = gexngenes(ismember(gexngenes, ds.dgexGenes));
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

idxnode = ismember(netpgenes, ds.mutGenes);
mutpgenes = union(mutpgenes, netpgenes(idxnode), 'stable');
idxnode = ismember(netngenes, ds.mutGenes);
mutngenes = union(mutngenes, netngenes(idxnode), 'stable');
idxnode = ismember(netpgenes, ds.cnvGenes);
cnvpgenes = union(cnvpgenes, netpgenes(idxnode), 'stable');
idxnode = ismember(netngenes, ds.cnvGenes);
cnvngenes = union(cnvngenes, netngenes(idxnode), 'stable');
idxnode = ismember(netpgenes, ds.dgexGenes);
gexpgenes = union(gexpgenes, netpgenes(idxnode), 'stable');
idxnode = ismember(netngenes, ds.dgexGenes);
gexngenes = union(gexngenes, netngenes(idxnode), 'stable');

idxdrug = strcmpi(ds.allDrugs, drug);
[~, idxsort] = sort(ds.AAMat(idxdrug, :));
%idxrest = find(ds.AAMat(idxdrug, idxsort)<=thaa(1)); % | ds.IC50Mat(idxdrug, idxsort)>=thic50(1));
%idxresp = find(ds.AAMat(idxdrug, idxsort)>=thaa(2)); % | ds.IC50Mat(idxdrug, idxsort)<=thic50(2));

idxmutp = cell2mat(cellfun(@(x) find(strcmp(x, ds.mutGenes)), mutpgenes, 'UniformOutput', 0));
idxcnvp = cell2mat(cellfun(@(x) find(strcmp(x, ds.cnvGenes)), cnvpgenes, 'UniformOutput', 0));
idxgexp = cell2mat(cellfun(@(x) find(strcmp(x, ds.dgexGenes)), gexpgenes, 'UniformOutput', 0));
idxmutn = cell2mat(cellfun(@(x) find(strcmp(x, ds.mutGenes)), mutngenes, 'UniformOutput', 0));
idxcnvn = cell2mat(cellfun(@(x) find(strcmp(x, ds.cnvGenes)), cnvngenes, 'UniformOutput', 0));
idxgexn = cell2mat(cellfun(@(x) find(strcmp(x, ds.dgexGenes)), gexngenes, 'UniformOutput', 0));


pmat = ds.mutMat([idxmutp; idxmutn], idxsort);
pmat = [pmat; ds.cnvMat([idxcnvp; idxcnvn], idxsort)];
pmat = [pmat; ds.dgexMat([idxgexp; idxgexn], idxsort)];

rowlabels = strcat([mutpgenes; mutngenes], '^{MUT}');
rowlabels = [rowlabels; strcat([cnvpgenes; cnvngenes], '^{CNV}')];
rowlabels = [rowlabels; strcat([gexpgenes; gexngenes], '^{GEX}')];
profs  = ds.AAMat(idxdrug, idxsort);



