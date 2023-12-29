function [err, corr] = rmse(pdp, pdpsim,idx)
err = zeros(length(idx),1);
corr = zeros(length(idx),1);
for i = 1:length(idx)
    ppp = 10*log10(pdp(:,idx(i))/max(pdp(:, idx(i))));
    pdpsim = 10*log10(abs(pdpsim));
    m = min(length(ppp), length(pdpsim));
    err(i) = sqrt(immse(ppp(1:m), pdpsim(1:m)));
    co = corrcoef(ppp(1:m), pdpsim(1:m));
    corr(i) = co(2,1);
end
end