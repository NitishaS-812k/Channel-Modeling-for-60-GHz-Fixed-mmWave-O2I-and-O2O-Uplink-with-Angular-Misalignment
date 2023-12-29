function [PDPm, PDPdB, idx, threshold] = preprocess_PDP(PDP)
PDP_norm = PDP/sum(PDP);
PDP_normdB = 10*log10(PDP_norm);
threshold = (-1)*mean(PDP_normdB);
idx = zeros(1, length(PDP));
maxm = max(PDP_normdB);
for i= 1:length(PDP)
    if PDP_normdB(i) < (maxm-threshold)
        PDP_normdB(i) = 0;
        idx(i) = i;
    end
end
PDP_norm(nonzeros(idx)) = 0;
PDPm = PDP_norm;
PDPdB = PDP_normdB;
idx = nonzeros(idx);
end