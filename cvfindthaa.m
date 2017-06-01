for alpha=0.9
     for maxn=50
        filename = sprintf('subnets/diffkfoldctrp_rng_def_%.2f_%d.csv', alpha, maxn);
        if(exist(filename, 'file')~=2)
            continue;
        end
        cvtab = readtable(filename);

        nofsacc = zeros(height(cvtab), 1);
        fsacc = nofsacc;
        diffacc = nofsacc;
        for i=1:height(cvtab)
            drug = cvtab.DRUG{i};

            tpnofs = cvtab{i, {sprintf('tpnofs')}};
            fpnofs = cvtab{i, {sprintf('fpnofs')}};
            tnnofs = cvtab{i, {sprintf('tnnofs')}};
            fnnofs = cvtab{i, {sprintf('fnnofs')}};
            nofsacc(i) = (tpnofs + tnnofs)/(tpnofs + tnnofs + fpnofs + fnnofs);

            tpfs = cvtab{i, {sprintf('tpfs')}};
            fpfs = cvtab{i, {sprintf('fpfs')}};
            tnfs = cvtab{i, {sprintf('tnfs')}};
            fnfs = cvtab{i, {sprintf('fnfs')}};
            fsacc(i) = (tpfs + tnfs)/(tpfs + tnfs + fpfs + fnfs);

            tpdiff = cvtab{i, {sprintf('tpdiff')}};
            fpdiff = cvtab{i, {sprintf('fpdiff')}};
            tndiff = cvtab{i, {sprintf('tndiff')}};
            fndiff = cvtab{i, {sprintf('fndiff')}};
            diffacc(i) = (tpdiff + tndiff)/(tpdiff + tndiff + fpdiff + fndiff);

        end
        idxnan = isnan(nofsacc) | isnan(diffacc) | isnan(fsacc);
        nofsacc = nofsacc(~idxnan);
        diffacc = diffacc(~idxnan);
        fsacc = fsacc(~idxnan);
        [h, p] = ttest(nofsacc, diffacc);
     fprintf('alpha = %.2f, maxn = %d, mean(NOFS) = %.2f, mean(DIFF) = %.2f, ACDIFF < ACCNOFS = %d, ACCDIFF > ACCNOFS = %d, p = %.2f\n', ...
         alpha, maxn, mean(nofsacc), mean(diffacc), sum(diffacc < nofsacc), sum(diffacc > nofsacc), p);
     end
     
end