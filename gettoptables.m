for i=1:length(cclenet.allDrugs)
    ccle_tbl{i} = topGenes(ccle_dgnet{i}, cclenet, net, 0.05, 1e-3);
    writetable(ccle_tbl{i}, sprintf('lists/ccle2_%s.csv', cclenet.allDrugs{i}));
end

for i=228:length(gdscnet.allDrugs)
    gdsc_tbl{i} = topGenes(gdsc_dgnet{i}, gdscnet, net, 0.05, 1e-5);
    writetable(gdsc_tbl{i}, sprintf('lists/allgdsc_%s.csv', gdscnet.allDrugs{i}));
end

for i=1:length(combnet.allDrugs)
    comb_tbl{i} = topGenes(comb_dgnet{i}, combnet, net, 1e-3, 1e-6);
    writetable(comb_tbl{i}, sprintf('lists/comb2_%s.csv', combnet.allDrugs{i}));
end
