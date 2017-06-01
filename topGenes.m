function toptbl = topGenes(dgnet, ds, net, thpv1, thpv2)
% Analysing kernel results and obtain a list of top selected features

    startmutp = dgnet.globpos.startMut;
    endmutp = dgnet.globpos.endMut;
    startcnvp = dgnet.globpos.startCNV;
    endcnvp = dgnet.globpos.endCNV;
    startgexp = dgnet.globpos.startGEX;
    endgexp = dgnet.globpos.endGEX;
    
    startmutn = dgnet.globneg.startMut;
    endmutn = dgnet.globneg.endMut;
    startcnvn = dgnet.globneg.startCNV;
    endcnvn = dgnet.globneg.endCNV;
    
    startgexn = dgnet.globneg.startGEX;
    endgexn = dgnet.globneg.endGEX;
    
    if(length(dgnet.simpos) < dgnet.globpos.size)
        dgnet.simpos = zeros(1, dgnet.globpos.size);
    end
    
    if(length(dgnet.simneg) < dgnet.globneg.size)
        dgnet.simneg = zeros(1, dgnet.globneg.size);
    end
    
    dgnet.simpos(isnan(dgnet.simpos)) = 0;
    dgnet.simneg(isnan(dgnet.simneg)) = 0;
    
    simdrugpos = dgnet.simpos*10e5;
    simdrugneg = dgnet.simneg*10e5;
	
    smut = [simdrugpos(startmutp:endmutp)' simdrugneg(startmutn:endmutn)'];
    scnv = [simdrugpos(startcnvp:endcnvp)' simdrugneg(startcnvn:endcnvn)'];
    sgex = [simdrugpos(startgexp:endgexp)' simdrugneg(startgexn:endgexn)'];
    
    smut = quantilenorm(smut);
    scnv = quantilenorm(scnv);
    sgex = quantilenorm(sgex);
    
    smut = smut(:, 1) - smut(:, 2);
    scnv = scnv(:, 1) - scnv(:, 2);
    sgex = sgex(:, 1) - sgex(:, 2);

    simnetp = dgnet.simpos(1:net.size)*10e5;
    simnetn = dgnet.simneg(1:net.size)*10e5;
    
    simnet = quantilenorm([simnetp' simnetn']);
    simnet = simnet(:, 1) - simnet(:, 2);

    zmut = zscore(smut);
    pv = (1-normcdf(abs(zmut), 0, 1));
    topmutpidx = zmut>0 & pv<=thpv1;
    topmutnidx = zmut<0 & pv<=thpv1;
    topmutgenesp = ds.mutGenes(topmutpidx);
    topmutgenesn = ds.mutGenes(topmutnidx);
    
    [~, idx] = sort(smut(topmutpidx), 'descend');
    topmutgenesp = topmutgenesp(idx);
    
    [~, idx] = sort(smut(topmutnidx));
    topmutgenesn = topmutgenesn(idx);
    
    zcnv = zscore(scnv);
    pv = (1-normcdf(abs(zcnv), 0, 1));
    toppidx = zcnv>0 & pv<=thpv1;
    topnidx = zcnv<0 & pv<=thpv1;
    topcnvgenesp = ds.cnvGenes(toppidx);
    topcnvgenesn = ds.cnvGenes(topnidx);
    
    [~, idx] = sort(scnv(toppidx), 'descend');
    topcnvgenesp = topcnvgenesp(idx);
    
    [~, idx] = sort(scnv(topnidx));
    topcnvgenesn = topcnvgenesn(idx);
    
    
    zgex = zscore(sgex);
    pv = (1-normcdf(abs(zgex), 0, 1));
    toppidx = zgex>0 & pv<=thpv1;
    topnidx = zgex<0 & pv<=thpv1;
    topgexgenesp = ds.dgexGenes(toppidx);
    topgexgenesn = ds.dgexGenes(topnidx);
    [~, idx] = sort(sgex(toppidx), 'descend');
    topgexgenesp = topgexgenesp(idx);
    
    [~, idx] = sort(sgex(topnidx));
    topgexgenesn = topgexgenesn(idx);
    
    
    znet = zscore(simnet);
    pv = (1-normcdf(abs(znet), 0, 1));
    toppidx = znet>0 & pv<=thpv2;
    topnidx = znet<0 & pv<=thpv2;
    
    topnetgenesp = net.nodes(toppidx);
    topnetgenesn = net.nodes(topnidx);
    
    [~, idx] = sort(simnet(toppidx), 'descend');
    topnetgenesp = topnetgenesp(idx);
    
    [~, idx] = sort(simnet(topnidx));
    topnetgenesn = topnetgenesn(idx);
    
    maxl = max([length(topmutgenesp); length(topmutgenesn); length(topcnvgenesp); length(topcnvgenesn); ...
        length(topgexgenesp); length(topgexgenesn); length(topnetgenesp); length(topnetgenesn)]);
    
    topmutgenesp{maxl+1} = '';
    topmutgenesn{maxl+1} = '';
    topcnvgenesp{maxl+1} = '';
    topcnvgenesn{maxl+1} = '';
         
    topgexgenesp{maxl+1} = '';
    topgexgenesn{maxl+1} = '';
    topnetgenesp{maxl+1} = '';
    topnetgenesn{maxl+1} = '';
    
    topmutgenesp = reshape(topmutgenesp, maxl+1, 1);
    topmutgenesn = reshape(topmutgenesn, maxl+1, 1);
    topcnvgenesp = reshape(topcnvgenesp, maxl+1, 1);
    topcnvgenesn = reshape(topcnvgenesn, maxl+1, 1);
    topgexgenesp = reshape(topgexgenesp, maxl+1, 1);
    topgexgenesn = reshape(topgexgenesn, maxl+1, 1);
    
    topnetgenesp = reshape(topnetgenesp, maxl+1, 1);
    topnetgenesn = reshape(topnetgenesn, maxl+1, 1);
    
    toptbl = table(topmutgenesp, topmutgenesn, topcnvgenesp, topcnvgenesn, topgexgenesp, topgexgenesn, topnetgenesp, topnetgenesn, ... '
        'VariableNames', {'MUTP', 'MUTN', 'CNVP', 'CNVN', 'GEXP', 'GEXN', 'NETP', 'NETN'});
    
end