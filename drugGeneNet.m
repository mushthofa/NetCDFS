function drugnet = drugGeneNet(ds, idxdrug, net, idxpos, idxneg)
    % Perform network integration and kernel calculation for responsive and
    % resistant sampels of a certain drug

% 
%     idxdrug = strcmpi(ds.allDrugs, drug);
%     if(strcmpi(criteria, 'aa'))
%         prof = ds.AAMat(idxdrug, :);
%     elseif(strcmpi(criteria, 'ic50'))
%         prof = ds.IC50Mat(idxdrug, :);
%     else
%         error('Wrong type of criteria. Must be "aa" or "ic50"!');
%     end
% 
%     idxneg = prof<=params(1);
%     idxpos = prof>=params(2);
%     drugnet.poscells = ds.cellNames(idxpos);
%     posprof = prof(idxpos);
%     drugnet.negcells = ds.cellNames(idxneg);
%     negprof = prof(idxneg);

    posprof = ds.IC50Mat(idxdrug, idxpos);
    negprof = ds.IC50Mat(idxdrug, idxneg);
    

    % Normalize weights;
%     if(strcmp(criteria, 'aa'))
        posprof = (posprof-min(posprof))./(max(posprof)-min(posprof));
        negprof = (negprof-min(negprof))./(max(negprof)-min(negprof));
        posprof = 1 - posprof;
%     else
%         posprof = (posprof-min(posprof))./(max(posprof)-min(posprof));
%         negprof = (negprof-min(negprof))./(max(negprof)-min(negprof));
%         posprof = 1 - posprof;
%     end

   
    drugnet.posprof = posprof;
    drugnet.negprof = negprof;

    % Build positive global network
    disp('Computing positive kernel...');
    globpos = buildNetMat(drugnet.posprof', ds.mutMat(:, idxpos), ds.mutGenes, ds.cnvMat(:, idxpos)~=0, ds.cnvGenes, ds.dgexMat(:, idxpos)~=0, ds.dgexGenes, net);
    kernelmat = kernel(globpos.mat, 'LEXP', 0.01);
    kernelmat = kernelNorm(kernelmat);
    simpos = mean(kernelmat(globpos.startSample:globpos.endSample, 1:end));
    globpos = rmfield(globpos, 'mat');



    disp('Computing negative kernel...');
    globneg = buildNetMat(drugnet.negprof', ds.mutMat(:, idxneg), ds.mutGenes, ds.cnvMat(:, idxneg)~=0, ds.cnvGenes, ds.dgexMat(:, idxneg)~=0, ds.dgexGenes, net);
    kernelmat = kernel(globneg.mat, 'LEXP', 0.01);
    kernelmat = kernelNorm(kernelmat);
    simneg = mean(kernelmat(globneg.startSample:globneg.endSample, 1:end));
    globneg = rmfield(globneg, 'mat');

    drugnet.simpos = simpos;
    drugnet.simneg = simneg;
    drugnet.globpos = globpos;
    drugnet.globneg = globneg;
end







