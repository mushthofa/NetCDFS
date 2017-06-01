function ds = restrictNet(ds, net)
% Restrict dataset ds on the network net
% i.e., eliminate all genes not on the network

    idxnav = ~ismember(ds.mutGenes, net.nodes);
    ds.mutGenes(idxnav) = [];
    ds.mutMat(idxnav, :) = [];

    idxnav = ~ismember(ds.cnvGenes, net.nodes);
    ds.cnvGenes(idxnav) = [];
    ds.cnvMat(idxnav, :) = [];

    idxnav = ~ismember(ds.gexGenes, net.nodes);
    ds.gexGenes(idxnav) = [];
    ds.gexMat(idxnav, :) = [];

%     idxnav = ~ismember(ds.dgexGenes, net.nodes);
%     ds.dgexGenes(idxnav) = [];
%     ds.dgexMat(idxnav, :) = [];

end

