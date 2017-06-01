pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
                         p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
pStart = .5;
lb = [0 -Inf -Inf 0 0];
ub = [1 Inf Inf Inf Inf];

fout = fopen('cvctrpthaa.csv', 'w');
fprintf(fout, 'DRUG,P,MU1,MU2,THLO,THHI,NUM0,NUM1,NUMREST\n');
for i=1:length(drugcvctrp)
    drug = drugcvctrp{i};
    idxd = strcmpi(ctrpnet.allDrugs, drug); 
    aa = ctrpnet.AAMat(idxd, :);
    aa(isnan(aa)) = [];
    
    muStart = quantile(aa,[.25 .75]);
    sigmaStart = sqrt(var(aa) - .25*diff(muStart).^2);
    start = [pStart muStart sigmaStart sigmaStart];
    options = statset('MaxIter', 1000, 'MaxFunEvals', 2000);
    paramEsts = mle(aa, 'pdf',pdf_normmixture, 'start',start, ...
                          'lower',lb, 'upper',ub, 'options',options);

    if(paramEsts(1) <= 0.05)
        thaa_lo = paramEsts(3);% - paramEsts(5);
        thaa_hi = paramEsts(3);% + paramEsts(5);
    elseif(paramEsts(1) >= 0.95)
        thaa_lo = paramEsts(2);% - paramEsts(4);
        thaa_hi = paramEsts(2);% + paramEsts(4);
    else
        thaa_lo = (paramEsts(1)*paramEsts(3) + (1-paramEsts(1))*paramEsts(2));
        thaa_hi = thaa_lo;
%         if(paramEsts(2) < paramEsts(3))
%             thaa_lo = paramEsts(2) ;
%             thaa_hi = paramEsts(3) ;
%          else
%             thaa_lo = paramEsts(3) ;
%             thaa_hi = paramEsts(2) ;
%         end
    end
    fprintf(fout, '%s,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%d,%d\n',drug, paramEsts(1),paramEsts(2), paramEsts(3), thaa_lo, thaa_hi, sum(aa<=thaa_lo), sum(aa>=thaa_hi), sum(aa<thaa_hi & aa>thaa_lo));
end
fclose(fout);