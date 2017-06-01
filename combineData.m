function ds = combineData(ds1, ds2)
% Combining two datasets into one
% e.g., GDSC + CCLE

    cl1 = ds1.cellNames;
    cl2 = ds2.cellNames;
    [cl, idx1, idx2] = union(cl1, cl2);
    tiss = [ds1.cellTissues(idx1); ds2.cellTissues(idx2)];
    idxcl1 = cellfun(@(x) find(strcmp(x, cl)), cl1);
    idxcl2 = cellfun(@(x) find(strcmp(x, cl)), cl2);

    % Combine mutation, take the maximum (OR)
    gmut1 = ds1.mutGenes;
    gmut2 = ds2.mutGenes;
    gmut = union(gmut1, gmut2);

    mutmat = zeros(length(gmut), length(cl));

    idxg1 = cellfun(@(x) find(strcmp(x, gmut)), gmut1);
    idxg2 = cellfun(@(x) find(strcmp(x, gmut)), gmut2);
    mutmat(idxg1, idxcl1) = max(mutmat(idxg1, idxcl1), ds1.mutMat);
    mutmat(idxg2, idxcl2) = max(mutmat(idxg2, idxcl2), ds2.mutMat);

    % Combine CNV, take the extreme
    gcnv1 = ds1.cnvGenes;
    gcnv2 = ds2.cnvGenes;
    gcnv = union(gcnv1, gcnv2);
    cnvmat = zeros(length(gcnv), length(cl));

    idxg1 = cellfun(@(x) find(strcmp(x, gcnv)), gcnv1);
    idxg2 = cellfun(@(x) find(strcmp(x, gcnv)), gcnv2);
    cnvmat(idxg1, idxcl1) = cnvmat(idxg1, idxcl1) + ds1.cnvMat;
    cnvmat(idxg2, idxcl2) = cnvmat(idxg2, idxcl2) + ds2.cnvMat;
    cnvmat = double(cnvmat>0) - double(cnvmat<0);

    % Combine GEX, intersect the gens, take the average,
    gex1 = ds1.gexGenes;
    gex2 = ds2.gexGenes;
    gex = intersect(gex1, gex2);
    gexmat = zeros(length(gex), length(cl));
    divmat = ones(length(gex), length(cl));

    idxg1 = cellfun(@(x) find(strcmp(x, gex1)), gex);
    idxg2 = cellfun(@(x) find(strcmp(x, gex2)), gex);
    gexmat(:, idxcl1) = gexmat(:, idxcl1) + ds1.gexMat(idxg1, :);
    gexmat(:, idxcl2) = gexmat(:, idxcl2) + ds2.gexMat(idxg2, :);
    divmat(:, idxcl1) = divmat(:, idxcl1) + 1;
    divmat(:, idxcl2) = divmat(:, idxcl2) + 1;

    gexmat = gexmat./divmat;

    % Save molecular data

    ds.mutMat = mutmat;
    ds.mutGenes = gmut;
    ds.cnvMat = cnvmat;
    ds.cnvGenes = gcnv;
    ds.gexMat = gexmat;
    ds.gexGenes = gex;
    ds.cellNames = cl;
    ds.cellTissues = tiss;

    % Combine profile data
    dr1 = ds1.allDrugs;
    dr2 = ds2.allDrugs;
    dr = union(dr1, dr2);

    idxdr1 = cellfun(@(x) find(strcmp(x, dr)), dr1);
    idxdr2 = cellfun(@(x) find(strcmp(x, dr)), dr2);
    idxcl1 = cellfun(@(x) find(strcmp(x, cl)), cl1);
    idxcl2 = cellfun(@(x) find(strcmp(x, cl)), cl2);

    ic1 = ds1.IC50Mat;
    ic2 = ds2.IC50Mat;
    icmat = nan(length(dr), length(cl));
    aa1 = ds1.AAMat;
    aa2 = ds2.AAMat;
    aamat = nan(length(dr), length(cl));

    % For IC50, we take the minimum
    icmat(idxdr1, idxcl1) = min(icmat(idxdr1, idxcl1), ds1.IC50Mat);
    icmat(idxdr2, idxcl2) = min(icmat(idxdr2, idxcl2), ds2.IC50Mat);

    % For AA, we take the maximum
    aamat(idxdr1, idxcl1) = max(aamat(idxdr1, idxcl1), ds1.AAMat);
    aamat(idxdr2, idxcl2) = max(aamat(idxdr2, idxcl2), ds2.AAMat);

    ds.IC50Mat = icmat;
    ds.AAMat = aamat;
    ds.allDrugs = dr;
end