function ds = marginalExp(ds, thdiff, thgexon, thstdgex)
% Compute marginal gene expression values
% 1, 0 or -1 depending on the deviations from tissue

    alltissues = unique(ds.cellTissues);
    dgexmat = zeros(size(ds.gexMat));
    for i=1:length(alltissues)
        idxtiss = ismember(ds.cellTissues, alltissues{i});
        gexmattis = ds.gexMat(:, idxtiss);

        gexmattop = repmat(quantile(gexmattis', thdiff)', 1, size(gexmattis, 2));
        gexmatbot = repmat(quantile(gexmattis', 1-thdiff)', 1, size(gexmattis, 2));
        gexmatavg = repmat(mean(gexmattis, 2), 1, size(gexmattis, 2));
        gexmatstd = repmat(std(gexmattis')', 1, size(gexmattis, 2));
        % Put 1 on genes that are >= thdiff quantile within tissue average and
        % whose tissue average >= thgexon
        % and vice versa for -1
        numd = sum(gexmattis<=gexmatbot, 2);
        numu = sum(gexmattis>=gexmattop, 2);
        selu = repmat(numu>numd | (numu.*numd == 1), 1, size(gexmattis, 2));
        seld = repmat(numd>numu | (numu.*numd == 1), 1, size(gexmattis, 2));
        
        dgexmattis = double(gexmatavg>=thgexon) .* double(gexmatstd>=thstdgex) .* (selu.*double(gexmattis>=gexmattop) - seld.*double(gexmattis<=gexmatbot));
        dgexmat(:, idxtiss) = dgexmattis;
    end
    
    idxnav = sum(dgexmat~=0, 2)==0;
    dgexmat(idxnav, :) = [];
    
    ds.dgexMat = dgexmat;
    ds.dgexGenes = ds.gexGenes(~idxnav);
    
end
