% Tissue-label matching between CCLE and GDSC
% Since we have shared cell lines between the two dataset
% we need to map tissue labels between the two, using the shared cell line
% and then propagate to the rest of the cell lines
% CCLE seems to be less specific about tissue label (e.g. does not
% differentiate between lung SCLC or NSCLC) while GDSC in general is more
% specific. Thus, to make a consistent mapping, we map each tissue label
% from GDSC to CCLE (with some manual changes below)

gdsc_cells = readtable('data/GDSC_Cells.csv');
ccle_cells = readtable('data/CCLE_Cells.csv');
% Standardize cell names
gdsccells = cellfun(@stdCellName, gdsc_cells.CELL_NAME, 'UniformOutput', 0);
cclecells = cellfun(@stdCellName, ccle_cells.Cell_Name, 'UniformOutput', 0);

% Find shared cell lines
gdsc_cells.CELL_NAME = gdsccells;
ccle_cells.Cell_Name = cclecells;
intcells = intersect(gdsccells, cclecells);
idxint1 = ismember(ccle_cells.Cell_Name, intcells);
idxint2 = ismember(gdsc_cells.CELL_NAME, intcells);

cellshared = table();
cellshared.CCLE_NAME = ccle_cells.Cell_Name(idxint1);
cellshared.CCLE_TISSUE = ccle_cells.site(idxint1);
cellshared = sortrows(cellshared);
gdscshared = table();
gdscshared.GDSC_NAME = gdsc_cells.CELL_NAME(idxint2);
gdscshared.GDSC_TISSUE = gdsc_cells.TISSUE2(idxint2);
gdscshared.CELL_ID = gdsc_cells.CELL_ID(idxint2);
gdscshared = sortrows(gdscshared);
cellshared.GDSC_NAME = gdscshared.GDSC_NAME;
cellshared.GDSC_TISSUE = gdscshared.GDSC_TISSUE;
cellshared.CELL_ID = gdscshared.CELL_ID;
clear gdscshared;

% Manually changed
cellshared.CCLE_TISSUE{16} = 'upper_aerodigestive_tract'; % A253: Previously salivary_gland, changed to u_a_t to make it consistent
cellshared.CCLE_TISSUE{23} = 'soft_tissue'; % A673: Inconsistent between CCLE and GDSC, we follow GDSC

% Mapping

tissuemap = containers.Map('KeyType', 'char','ValueType', 'char');
for i=1:height(cellshared)
     ccle_t = lower(cellshared.CCLE_TISSUE{i});
     gdsc_t = lower(cellshared.GDSC_TISSUE{i});
     if(~strcmp(gdsc_t, ccle_t))
         if(ismember(gdsc_t, keys(tissuemap)) && ~strcmp(tissuemap(gdsc_t), ccle_t))
             % This should not occur (i.e., after manual changes, we don't
             % have any inconsistent mapping
             fprintf('Line %d: %s was previously mapped to %s, now mapped to %s\n',i, ccle_t, tissuemap(ccle_t), ccle_t);
         end
         if(~ismember(gdsc_t, keys(tissuemap)))
             tissuemap(gdsc_t) = ccle_t;
         end
     end
end

idxgdscrep = ismember(lower(gdsc_cells.TISSUE2), keys(tissuemap));
gdsc_tissue_rep = gdsc_cells.TISSUE2(idxgdscrep);
gdsc_cells.MAPPEDTISS = lower(gdsc_cells.TISSUE2);
gdsc_cells.MAPPEDTISS(idxgdscrep) = values(tissuemap, lower(gdsc_cells.TISSUE2(idxgdscrep)));

