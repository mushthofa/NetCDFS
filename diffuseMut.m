function ds = diffuseMut(ds, net, alpha, minscore)
% Diffuse mutation signals in dataset ds accros the network net
% 

    maxn = 50;
    eps = 1e-5;
    idxmut = cell2mat(values(net.name2node, ds.mutGenes));

    nSample = size(ds.mutMat, 2);
    diffmutMat = zeros(net.size, nSample);
    
    diffmutMat(idxmut, :) = ds.mutMat;
    
    diffmutMat = netSmooth(diffmutMat', net.mat, alpha, eps, maxn)';
    
    % Only takes the genes that have at least one sample with value >=
    % minscore
    dm = sum(diffmutMat' >= minscore);
    idxdm = dm > 0;
    diffmutMat = diffmutMat(idxdm, :);
    
    ds.mutMat = double(diffmutMat>minscore);
    ds.mutGenes = net.nodes(idxdm);
    
end