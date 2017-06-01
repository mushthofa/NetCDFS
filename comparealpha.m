for alpha=0.9
    for maxn = 30
        for ntrees = 30
            filename = sprintf('./subnets/diffkfold_rng_def_%.2f_%d_%d.csv', alpha, maxn, ntrees);
            %disp(filename);
            if(exist(filename, 'file') ~= 2)
                continue;
            end
                
            cvtab = readtable(filename);
            if(height(cvtab) < 72)
                continue;
            end
            nofsacc = (cvtab.TPNOFS + cvtab.TNNOFS)./(cvtab.TPNOFS + cvtab.TNNOFS + cvtab.FPNOFS + cvtab.FNNOFS);
            fsacc = (cvtab.TPFS + cvtab.TNFS)./(cvtab.TPFS + cvtab.TNFS + cvtab.FPFS + cvtab.FNFS);
            diffacc = (cvtab.TPDIFF + cvtab.TNDIFF)./(cvtab.TPDIFF + cvtab.TNDIFF + cvtab.FPDIFF + cvtab.FNDIFF);
  
            nofsposrec = (cvtab.TPNOFS)./(cvtab.TPNOFS + cvtab.FNNOFS);
            fsposrec = (cvtab.TPFS)./(cvtab.TPFS + cvtab.FNFS);
            diffposrec = (cvtab.TPDIFF)./(cvtab.TPDIFF + cvtab.FNDIFF);

            [h, p] = ttest(diffacc, nofsacc);
            fprintf('alpha = %.2f, maxn = %d, ntrees = %d, DIFF > NOFS = %d, DIFF == NOFS = %d, DIFF < NOFS = %d, ttest pval = %f\n', alpha, maxn, ntrees, sum(diffacc > nofsacc), sum(diffacc==nofsacc), sum(diffacc<nofsacc), p);

        %     [h, p] = ttest(diffposrec, nofsposrec);
        %     fprintf('alpha = %.2f, DIFF > NOFS = %d, DIFF == NOFS = %d, ttest pval = %f\n', alpha, sum(diffposrec > nofsposrec), sum(diffposrec==nofsposrec), p);
        end
    end
end
