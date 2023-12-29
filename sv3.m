function [g, G, GoF,l, L] = sv3(PDP, tau)
%g: [1, num_clusters] shaped array with intra cluster slopes
%G: slope of inter cluster line
%GoF: Goodness of fit for intra cluster slopes

%constants
num_clusters = 3;
g = zeros(1, num_clusters);
l = zeros(1, num_clusters);
cluster_indices = ones(1, num_clusters);
GoF = zeros(1, num_clusters);
min_distance = 10;
ts = zeros(1, num_clusters);
%threshold PDP to obtain significant peaks
[pdp, pdpdb, idx, thresh] = preprocess_PDP(PDP);

%nonzero indices
idc = setdiff(1:length(pdpdb), idx);

%find indices of clusters
[~, ppos] = findpeaks(pdpdb(idc), 'MinPeakHeight', (-1)*thresh, 'MinPeakDistance', min_distance);

%in case three peaks are identified, pick the larger of the two
if length(ppos)< num_clusters
    cluster_indices(2) = ppos(1);
    cluster_indices(3) = ppos(2);
elseif length(ppos) >= num_clusters
    cluster_indices(2) = ppos(2);
    cluster_indices(3) = ppos(3);
end

%fit options
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );

%intra cluster slope
cluster_1 = pdp(1:cluster_indices(2)); %cluster 1
t1 = tau(1:cluster_indices(2));
idx1 = find(cluster_1 == 0);
idx1c = setdiff(1:length(t1), idx1);
t11 = t1(idx1c);
cluster_1n = cluster_1/max(cluster_1);
cluster_1dB = 10*log10(cluster_1n);
[fr1, gof1] = clusterslope(cluster_1dB, t1', ft, opts);
l(1) = 1/mean(diff(t11));
plot(fr1, t1, cluster_1)

%cluster 2
cluster_2 = pdp(cluster_indices(2):cluster_indices(3));
t2 = tau(cluster_indices(2):cluster_indices(3));
idx2 = find(cluster_2 == 0);
idx2c = setdiff(1:length(t2), idx2);
t22 = t2(idx2c);
cluster_2n = cluster_2/max(cluster_2);
cluster_2dB = 10*log10(cluster_2n);
[fr2, gof2] = clusterslope(cluster_2dB, t2', ft, opts);
l(2) = 1/mean(diff(t22));

%cluster 3
cluster_3 = pdp(cluster_indices(3):length(pdp));
t3 = tau(cluster_indices(3):length(pdp));
cluster_3n = cluster_3/max(cluster_3);
idx3 = find(cluster_3 == 0);
idx3c = setdiff(1:length(t3), idx3);
t33 = t3(idx3c);
cluster_3dB = 10*log10(cluster_3n);
[fr3, gof3] = clusterslope(cluster_3dB, t3', ft, opts);
l(3) = 1/mean(diff(t33));

m = max([pdp(1), pdp(cluster_indices(2)), pdp(cluster_indices(3))]);
pdp_c = pdp/m;
g(1) = fr1.p1; g(2) = fr2.p1; g(3) = fr3.p1;
GoF(1) = gof1.rmse; GoF(2) = gof2.rmse; GoF(3) = gof3.rmse;

%inter cluster slope
x = 1:length(pdp);
x(1) = []; x(cluster_indices(2)) = []; x(cluster_indices(3)) = [];
GG = fit(tau', pdp_c , ft, 'Exclude', x);
G = GG.p1;
ts(1) = tau(cluster_indices(1)); ts(2) = tau(cluster_indices(2)); ts(3) = tau(cluster_indices(3));
L = 1/mean(diff(ts));
end