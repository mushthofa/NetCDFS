function ds = loadCTRP(entrezmap)
%  LOADCTRP Summary of this function goes here
%  Function to load CTRP data

    ccle_mut = parse_gct('data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct');
    ccle_gex = parse_gct('data/CCLE_Expression_Entrez_2012-09-29.gct', 'checkid', false);
   
    ccle_cells = readtable('data/CCLE_Cells.csv');
    cclecellnames = cellfun(@stdCellName, ccle_cells.Cell_Name, 'UniformOutput', 0);
    ccle_cells.Cell_Name = cclecellnames;
    
    % To make it consistent with GDSC, we change the tissue-label of these
    % cell lines manually, see gdsc_ccle_tissuemap.m
    ccle_cells.site{23} = 'upper_aerodigestive_tract'; %A253: Previously salivary_gland, changed to u_a_t to make it consistent
    ccle_cells.site{30} = 'soft_tissue'; % A673: Inconsistent between CCLE and GDSC, we follow GDSC
    
    % Match CCLEName in the mutation matrix with the reference table
    
    idxcells = cellfun(@(x) find(strcmp(x, ccle_cells.CCLEName)), ccle_mut.cid, 'UniformOutput', 0);
    idxnotfound = find(cellfun(@isempty, idxcells));
    
    % For some reasons, the cell line names that start with numbers get appended
    % with X in front in the mutation labels, we need to remove this
    
    for i=1:length(idxnotfound)
        if(ccle_mut.cid{idxnotfound(i)}(1) == 'X')
            ccle_mut.cid{idxnotfound(i)} = ccle_mut.cid{idxnotfound(i)}(2:end);
        end
    end
    
    
    % Match again
    idxcells = cellfun(@(x) find(strcmp(x, ccle_cells.CCLEName)), ccle_mut.cid, 'UniformOutput', 0);
    idxnotfound = find(cellfun(@isempty, idxcells));
    
    % Remove samples that are not matched to the Cell lines table
    
    ccle_mut.cid(idxnotfound) = [];
    ccle_mut.mat(:, idxnotfound) = [];
    
    % Do the same for GEX matrix
    
    idxcells = cellfun(@(x) find(strcmp(x, ccle_cells.CCLEName)), ccle_gex.cid, 'UniformOutput', 0);
    idxnotfound = find(cellfun(@isempty, idxcells)); % This results in an empty result
    
    ccle_gex.cid(idxnotfound) = [];
    ccle_gex.mat(:, idxnotfound) = [];
    
    
    % Find intersections of cell lines from MUT and GEX
    icl = intersect(ccle_mut.cid, ccle_gex.cid);
    
    idxc1 = ismember(ccle_mut.cid, icl);
    idxc2 = ismember(ccle_gex.cid, icl);
    
    % Consider only cell lines in the intersection
    ccle_mut.cid(~idxc1) = [];
    ccle_mut.mat(:, ~idxc1) = [];
    
    ccle_gex.cid(~idxc2) = [];
    ccle_gex.mat(:, ~idxc2) = [];
    

    
    % For some reasons, in the CCLE GEX matrix, NCIH292_LUNG appears twice
    idx2 = find(strcmp(ccle_gex.cid, 'NCIH292_LUNG'));
    
    % Since their correlation is high ~0.98, we just take the average of
    % the two
    
    ccle_gex.mat(:, idx2(1)) = mean(ccle_gex.mat(:, idx2), 2);
    ccle_gex.mat(:, idx2(2)) = [];
    ccle_gex.cid(idx2(2)) = [];

    % Reorder samples in MUT and GEX matrices
    idxc1 = cellfun(@(x) find(strcmp(x, ccle_mut.cid)), icl);
    idxc2 = cellfun(@(x) find(strcmp(x, ccle_gex.cid)), icl);
    
    ccle_mut.cid = ccle_mut.cid(idxc1);
    ccle_mut.mat = ccle_mut.mat(:, idxc1);
    ccle_gex.cid = ccle_gex.cid(idxc2);
    ccle_gex.mat = ccle_gex.mat(:, idxc2);
    
    % Add standardissed cell names and tissues to the MUT and GEX data
    idxcells = cell2mat(cellfun(@(x) find(strcmp(x, ccle_cells.CCLEName)), ccle_mut.cid, 'UniformOutput', 0));
    ccle_mut.CELL_NAME = ccle_cells.Cell_Name(idxcells);
    ccle_mut.TISSUE = ccle_cells.site(idxcells);
    


   
    % Delete tissues with low samples
    minsample = 10;
    tt=tabulate(ccle_mut.TISSUE);
    tisslow = tt(cell2mat(tt(:, 2))<minsample, 1);
    
    idxlow = ismember(ccle_mut.TISSUE, tisslow);
    ccle_mut.cid(idxlow) = [];
    ccle_mut.mat(:, idxlow) = [];
    ccle_mut.CELL_NAME(idxlow) = [];
    ccle_mut.TISSUE(idxlow) = [];
    
    ccle_gex.cid(idxlow) = [];
    ccle_gex.mat(:, idxlow) = [];
    ccle_gex.CELL_NAME = ccle_mut.CELL_NAME;
    ccle_gex.TISSUE = ccle_mut.TISSUE;
   
    % Process mutation data
    idxmut = ~cellfun(@isempty, strfind(ccle_mut.rid, '_MUT'));
    idxdel = ~cellfun(@isempty, strfind(ccle_mut.rid, '_DEL'));
    idxamp = ~cellfun(@isempty, strfind(ccle_mut.rid, '_AMP'));

    mutmat = ccle_mut.mat(idxmut, :);
    delmat = ccle_mut.mat(idxdel, :);
    ampmat = ccle_mut.mat(idxamp, :);

    allnames = cellfun(@(x) x(1:end-4), ccle_mut.rid, 'UniformOutput', false);
    namemut = allnames(idxmut);
   
    nameamp = allnames(idxamp);
    
    % Mutation data has special 'manually curated' mutation set for several
    % genes (e.g., KRAS G12-13, BRAF V600E) in addition to the 'normal'
    % mutation set for those genes
    % We opted to integrate this two types of mutation into one for each
    % gene
    
    idxmanmut = find(cellfun(@(x) ~isempty(strfind(x, '.')), namemut));
    for i=1:length(idxmanmut)
        ss = strsplit(namemut{idxmanmut(i)}, '.');
        genename = ss{1};
        idxg = strcmp(namemut, genename);
        mutmat(idxg, :) = mutmat(idxg, :) | mutmat(idxmanmut(i), :);
    end
    mutmat(idxmanmut, :) = [];
    namemut(idxmanmut) = [];
    
    gexeid = cellfun(@gexid2eid, ccle_gex.rid);
    
    idxav = ismember(gexeid, cell2mat(keys(entrezmap.entrez2name)));
    namegex = values(entrezmap.entrez2name, num2cell(gexeid(idxav)));
    ccle_gex.mat = ccle_gex.mat(idxav, :);
    
    
    % Store all cell line molecular data
    ds.mutMat = mutmat;
    ds.mutGenes = namemut;
    ds.cnvMat = ampmat - delmat;
    ds.cnvGenes = nameamp;
    ds.gexMat = ccle_gex.mat;
    ds.gexGenes = namegex;
    ds.cellNames = ccle_mut.CELL_NAME;
    ds.cellTissues = ccle_mut.TISSUE;
    
    %Process the drug profiles;
    
    
    disp('Reading CTRP drug profiles...');
    ctrp_drugs = readtable('data/v20.meta.per_compound.txt', 'Delimiter', '\t');
    ctrp_cells = readtable('data/v20.meta.per_cell_line.txt', 'Delimiter', '\t');
    ctrp_exp = readtable('data/v20.meta.per_experiment.txt', 'Delimiter', '\t');
    ctrp_auc = readtable('data/v20.data.curves_post_qc.txt', 'Delimiter', '\t');
    
  
   
    
    % Map between ids (experiment, ccl, cmpound)
    
    cclid2cellname = containers.Map('KeyType', 'double', 'ValueType', 'char');
    for i=1:height(ctrp_cells)
        cclid2cellname(ctrp_cells.master_ccl_id(i)) = ctrp_cells.ccl_name{i};
    end
    
    exp2cellname = containers.Map('KeyType', 'double', 'ValueType', 'char');
    for i=1:height(ctrp_exp)
        ccl_id = ctrp_exp.master_ccl_id(i);
        if(ismember(ccl_id, cell2mat(keys(cclid2cellname))))
            cellname = cclid2cellname(ccl_id);
            exp2cellname(ctrp_exp.experiment_id(i)) = cellname;
        end
    end
    
    ctrp_drugs.cpd_name = upper(ctrp_drugs.cpd_name);
    cpdid2name = containers.Map('KeyType', 'double', 'ValueType', 'char');
    for i=1:height(ctrp_drugs)
        cpdid2name(ctrp_drugs.master_cpd_id(i)) = ctrp_drugs.cpd_name{i};
    end
    
    % Only use experiments with known cell lines
    
    idxnav = ~ismember(ctrp_auc.experiment_id, cell2mat(keys(exp2cellname)));
    ctrp_auc(idxnav, :) = [];
    
    ctrp_auc.cellnames = values(exp2cellname, table2cell(ctrp_auc(:, 1)));
    
    % Only use experiments with CCLE cell lines
    
    idxnav = ~ismember(ctrp_auc.cellnames, ds.cellNames);
    ctrp_auc(idxnav, :) = [];
    
    
    % Only use experiments with mapped compound name
    idxnav = ~ismember(ctrp_auc.master_cpd_id, cell2mat(keys(cpdid2name)));
    ctrp_auc(idxnav, :) = [];
    
	ctrp_auc.cpd_name = values(cpdid2name, table2cell(ctrp_auc(:, end-1)));
    
    
    
    alldrugs = unique(ctrp_auc.cpd_name);
    avcells = unique(ctrp_auc.cellnames);
    
    
    % Delete cell lines data not used in the experiments
    idxnav = ~ismember(ds.cellNames, avcells);
    ds.mutMat(:, idxnav) = [];
    ds.cnvMat(:, idxnav) = [];
    ds.gexMat(:, idxnav) = [];
    ds.cellNames(idxnav) = [];
    ds.cellTissues(idxnav) = [];

    nDrugs = length(alldrugs);
    nCells = length(ds.cellNames);
    
    AAMat = nan(nDrugs, nCells);
    for i=1:nDrugs
        idxd = strcmp(ctrp_auc.cpd_name, alldrugs{i});
        idxc = cellfun(@(x) find(strcmp(ds.cellNames, x)), ctrp_auc.cellnames(idxd));
        AAMat(i, idxc) = ctrp_auc.area_under_curve(idxd);
    end
 
    ds.allDrugs = alldrugs;
    % Reverse AA from AUC, and map to 0-1
    ds.AAMat = (max(max(AAMat)) - AAMat)./max(max(AAMat));
    
    

end


