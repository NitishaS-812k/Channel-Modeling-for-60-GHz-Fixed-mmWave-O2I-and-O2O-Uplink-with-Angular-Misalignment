%calculating SV Parameters for O2I (dataset1) and O2O(dataset3)
%dataset 1
%npos = SizemagSdBV(2);
gamma1 = zeros(npos, 2);
GAMMA1 = zeros(npos, 1);
lambda1 = zeros(npos, 2);
LAMBDA1 = zeros(npos, 1);
GoF1 = zeros(npos, 2);
PADPp1 = PADP(timeeval,:);
CF11 = zeros(npos, 2);
CF21 = zeros(npos, 2);
for i = 1:npos
    [g1, G1, gof1, cf1, cf2, l1, L1] = sv1(PADPp1(:,i), tau);
    gamma1(i, :) = g1;
    GAMMA1(i) = G1;
    GoF1(i, :) = gof1;
    CF11(i, :) = cf1;
    CF21(i, :) = cf2;
    lambda1(i, :) = l1;
    LAMBDA1(i) = L1;
end

%dataset 2
%gamma2 = zeros(SizemagSdBV2(2), 2);
%GAMMA2 = zeros(SizemagSdBV2(2), 1);
%lambda2 = zeros(SizemagSdBV2(2), 2);
%LAMBDA2 = zeros(SizemagSdBV2(2), 1);
%GoF2 = zeros(SizemagSdBV2(2), 2);
%CF12 = zeros(SizemagSdBV2(2), 2);
%CF22 = zeros(SizemagSdBV2(2), 2);
%PADP_p2 = PADP2(timeeval2, :);
%for i = 1:SizemagSdBV2(2)
    %[g2, G2, gof2, cf1, cf2, l2, L2] = sv2(PADP_p2(:,i), timeee2);
    %gamma2(i, :) = g2;
    %GAMMA2(i) = G2;
    %GoF2(i, :) = gof2;
    %CF12(i, :) = cf1;
    %CF22(i, :) = cf2;
    %lambda2(i,:) = l2;
    %LAMBDA2(i) = L2;
%end

%dataset 3
gamma3 = zeros(SizemagSdBV3(2), 3);
GAMMA3 = zeros(SizemagSdBV3(2), 1);
GoF3 = zeros(SizemagSdBV3(2), 3);
lambda3 = zeros(SizemagSdBV3(2), 3);
LAMBDA3 = zeros(SizemagSdBV3(2), 1);
PADP_p3 = PADP3(timeeval3,:);
for i = 1:SizemagSdBV3(2)
    [g3, G3, gof3, l3, L3] = sv3(PADP_p3(:,i), timeee3);
    gamma3(i, :) = g3;
    GAMMA3(i) = G3;
    GoF3(i, :) = gof3;
    lambda3(i,:) = l3;
    LAMBDA3(i) = L3;
end

cl11 = gamma1(:, 1); %cluster 1
cl21 = gamma1(:, 2); %cluster 2
l11 = lambda1(:, 1);
l21 = lambda1(:, 2);
cl13 = gamma3(:, 1); %cluster 1
cl23 = gamma3(:, 2); %cluster 2
cl33 = gamma3(:, 3); %cluster 3
l13 = lambda3(:, 1);
l23 = lambda3(:, 2);
l33 = lambda3(:, 3);
