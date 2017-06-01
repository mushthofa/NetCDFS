intmut = intersect(gdsc.mutGenes, ctrp.mutGenes);
idxgdscint = cellfun(@(x)find(strcmp(x, gdsc.mutGenes)), intmut);
idxctrpint = cellfun(@(x)find(strcmp(x, ctrp.mutGenes)), intmut);

intcell = intersect(gdsc.cellNames, ctrp.cellNames);
idxgdscintc = cellfun(@(x)find(strcmp(x, gdsc.cellNames)), intcell);
idxctrpintc = cellfun(@(x)find(strcmp(x, ctrp.cellNames)), intcell);

mutmatg = gdsc.mutMat(idxgdscint, idxgdscintc);
mutmatc = ctrp.mutMat(idxctrpint, idxctrpintc);

HeatMap(mutmatg');
HeatMap(mutmatc');

missedg = mutmatg & ~mutmatc;
missedc = ~mutmatg & mutmatc;

nmg = sum(missedg, 2);
nmc = sum(missedc, 2);
figure; hist(nmg, 30);
figure; hist(nmc, 30);

