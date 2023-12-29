function [fr,gof] = clusterslope(cluster, time, ft, opts)
%cluster: the cluster array
%time: the time array
%ft: fittype
%opts: fitoptions

idxcl = find(cluster == -Inf);
weights = (abs(1./cluster)); 
t = find(weights == Inf);
weights(t) = 0.5;
opts.Exclude = idxcl;  %excluding these indices because they're zero
opts.Weights = weights;
[fr, gof] = fit(time, cluster ,ft, opts);
end
