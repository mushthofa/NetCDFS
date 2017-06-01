for alpha=0.5:0.05:0.95
    for maxn = 10:10:30
        for ntrees = 50:50:100
            filename = sprintf('./subnets/diffkfold_rng_def_%.2f_%d_%d.csv', alpha, maxn, ntrees);
            disp(filename);
            if(exist(filename, 'file') == 2)
                continue;
            end
            gensubnet;
        end
    end
end