function [g, G, GoF, confint1, confint2, l, L] = sv1(PDP, tau)
%g: [1, num_clusters] shaped array with intra cluster slopes
%G: slope of inter cluster line
%GoF: Goodness of fit for intra cluster slopes
%confint1: confidence intervals for slope of cluster 1
%confint2: confidence intervals of slope of cluster 2

%constants
num_clusters = 2;
g = zeros(1, num_clusters);
l = zeros(1, num_clusters);
GoF = zeros(1, num_clusters);
min_distance = 11;

%threshold PDP to obtain significant peaks
[pdp, pdpdb, idx, thresh] = preprocess_PDP(PDP);

%nonzero indices
idc = setdiff(1:length(pdpdb), idx);

%find indices of clusters
[pos, ppos] = findpeaks(pdpdb(idc), 'MinPeakHeight', (-1)*thresh, 'MinPeakDistance', min_distance);

%in case three peaks are identified, pick the larger of the two
if length(ppos) > num_clusters
    ix = find(pdpdb == max(pos(2:length(pos))));
else
    ix = ppos(2);
end

%fit options
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );

%intra cluster slope
cluster_1 = pdp(1:ix); %cluster 1
t1 = tau(1:ix);
idx1 = find(cluster_1 == 0);
idx1c = setdiff(1:length(t1), idx1);
t11 = t1(idx1c);
cluster_1n = cluster_1/max(cluster_1);
cluster_1dB = 10*log10(cluster_1n);
[fr1, gof1] = clusterslope(cluster_1dB, t1', ft, opts);
conf1 = confint(fr1);
confint1 = conf1(1, :);
l(1) = 1/mean(diff(t11));
plot(fr1, t1, cluster_1)

%cluster 2
cluster_2 = pdp(ix:length(pdpdb));
t2 = tau(ix:length(pdpdb));
idx2 = find(cluster_1 == 0);
idx2c = setdiff(1:length(t2), idx2);
t22 = t2(idx2c);
cluster_2n = cluster_2/max(cluster_2);
cluster_2dB = 10*log10(cluster_2n);
[fr2, gof2] = clusterslope(cluster_2dB, t2', ft, opts);
conf2 = confint(fr2);
confint2 = conf2(1,:);
l(2) = 1/mean(diff(t22));

m = max(pdp(1), pdp(ix));
pdp_c = pdp/m;
g(1) = fr1.p1; g(2) = fr2.p1;
GoF(1) = gof1.rmse; GoF(2) = gof2.rmse;

%inter cluster slope
x = 1:length(pdp);
x(1) = []; x(ix) = [];
GG = fit(tau', pdp_c , ft, 'Exclude', x);
G = GG.p1;
L = 1/tau(ix);
end