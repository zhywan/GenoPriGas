%Calculate the LR statistic for a single SNP
function lr=singlesnplr(genotype, poolfreq, reffreq)

switch genotype
    case 0
        lr = 2 * log((1-poolfreq)/(1-reffreq));
    case 1
        lr = log((1-poolfreq)/(1-reffreq)) + log (poolfreq/reffreq);
        %lr = log ((poolfreq * (1-poolfreq))/(reffreq * (1-reffreq)));%alternative way        
    case 2
        lr = 2 * log (poolfreq/reffreq);
    case -1
        lr = 0;
    otherwise
        lr = 0;
        disp('Unknown genotype');
end
