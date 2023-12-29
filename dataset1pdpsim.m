
%SV Parameters gamma(g), Gamma(G), lambda(l) and Lambda(L) sorted by
%elevation(0, 5, -5)
g1_0 = cl11(el0);     g2_0 = cl21(el0);   l1_0 = l11(el0);     l2_0 = l21(el0);
g1_n5 = cl11(eln5);  g2_n5 = cl21(eln5); l1_n5 = l11(eln5);    l2_n5 = l21(eln5);
g1_p5 = cl11(elp5);  g2_p5 = cl21(elp5); l1_p5 = l11(elp5);    l2_p5 = l21(elp5);

G1_0 = GAMMA1(el0);     L1_0 = LAMBDA1(el0);
G1_n5 = GAMMA1(eln5);  L1_n5 = LAMBDA1(eln5);
G1_p5 = GAMMA1(elp5);  L1_p5 = LAMBDA1(elp5);

%misalignment of 10
mis10_0 = find(A0 >= -10 & A0 <= 10);
los = mis10_0(3);
mis10_0(3) = [];
mis10_n5 = find(An5 >=- 10 & An5 <= 10); 
mis10_p5 = find(Ap5 >=- 10 & Ap5 <= 10); 

g1_m10 = 10^(mean([ g1_n5(mis10_n5)', g1_0(mis10_0)', g1_p5(mis10_p5)'])/10);
g2_m10 = 10^(mean([g2_n5(mis10_n5)', g2_0(mis10_0)', g2_p5(mis10_p5)'])/10);
%g2_m10 = 0.4;
l1_m10 = mean([l1_n5(mis10_n5)',l1_0(mis10_0)', l1_p5(mis10_p5)']);
l2_m10 = mean([l2_n5(mis10_n5)',l2_0(mis10_0)', l2_p5(mis10_p5)']);
G1_m10 = 10^(mean([G1_0(mis10_0)', G1_n5(mis10_n5)', G1_p5(mis10_p5)'])/10);
%G1_m10 = 0.75;
L1_m10 = mean([ L1_n5(mis10_n5)',L1_0(mis10_0)', L1_p5(mis10_p5)']);

%misalignment of 10-25 
mis25_0 = [find(A0 > 10 & A0 <= 25), find(A0 < -10 & A0 >= -25)];
mis25_n5 = [find(An5 > 10 & An5 <= 25), find(An5 < -10 & An5 >= -25)];
mis25_p5 = [find(Ap5 > 10 & Ap5 <= 25), find(Ap5 < -10 & Ap5 >= -25)];

g1_m25 = 10^(mean([g1_n5(mis25_n5)', g1_0(mis25_0)', g1_p5(mis25_p5)'])/10);
g2_m25 = 10^(mean([g2_n5(mis25_n5)',g2_0(mis25_0)', g2_p5(mis25_p5)'])/10);
%g2_m25 = 0.5;
l1_m25 = mean([l1_0(mis25_0)', l1_n5(mis25_n5)', l1_p5(mis25_p5)']);
l2_m25 = mean([l2_0(mis25_0)', l2_n5(mis25_n5)', l2_p5(mis25_p5)']);
G1_m25 = 10^(mean([G1_n5(mis25_n5)',G1_0(mis25_0)', G1_p5(mis25_p5)'])/10);
%G1_m25 = 0.8;
L1_m25 = mean([ L1_n5(mis25_n5)',L1_0(mis25_0)', L1_p5(mis25_p5)']);

%LOS
los01 = 10^(g1_0(los)/10); %gamma1 for los
%los01 = 0.21; %0.1
los02 = 10^(g2_0(los)/10); %gamma2 for los
%los02 = 0.58;%0.3;
l_los1 =  l1_0(los);%lambda 1 for los
l_los2 = l2_0(los);%lambda 2 for los
G_los = 10^(G1_0(los)/10);
%G_los = 0.45;%0.6;
L_los = L1_0(los);
[p00,t00,~,np00] = pdp_svsim_og(L_los,[l_los1, l_los2],G_los,[los01,los02],2,1,1,3,0,0);
figure(1)
plot(t00(1:np00(1),1), 10*log10(abs(p00(1:np00(1),1))))
hold on
p0 = 10*log10(PDP_0(:,6)/max(PDP_0(:, 6)));
plot(tau, p0)
hold off
legend('Simulated', '0_0')
title('LOS')

%test PDPs
%for PDPs in 0-10 misalignment range
[hs10,ts10,t0s10,np10] = pdp_svsim_og(L1_m10,[l1_m10, l2_m10],G1_m10,[g1_m10,g2_m10],2,1,1,3,0,0);
mean_del_10 = sum(abs(hs10).* ts10)/sum(abs(hs10));
rms_del_spr_10 = sqrt(sum(abs(hs10).*(ts10 - mean_del_10).^2)/sum(abs(hs10)));
figure(2)
plot(ts10(1:np10(1),1), 10*log10(abs(hs10(1:np10(1),1))))
hold on
p0_p5 = 10*log10(PDP_0(:,7)/max(PDP_0(:, 7)));
pn5_n5 = 10*log10(PDP_n5(:, 5)/max(PDP_n5(:,5)));
pp5_10 = 10*log10(PDP_p5(:, 8)/max(PDP_p5(:,8)));
plot(tau, p0_p5)
plot(tau, pn5_n5)
plot(tau, pp5_10)
hold off
legend('Simulated','0 0','-5 -5','+5 +10')
title('misalignment between 0-10 degrees')

%calculate rmse
[rmse_mis10_0, corr_mis10_0] = rmse(PDP_0, hs10, mis10_0);
[rmse_mis10_p5, corr_mis10_p5] = rmse(PDP_p5, hs10, mis10_p5);
[rmse_mis10_n5, corr_mis10_n5] = rmse(PDP_n5, hs10, mis10_n5);

%los rmse
%rmse_los = sqrt(immse(p0, p00));

%for PDPs in 10-25 misalignment range
[hs25,ts25,t0s25,np25] = pdp_svsim_og(L1_m25,[l1_m25, l2_m25],G1_m25,[g1_m25,g2_m25],2,1,1,3,0,0);
mean_del_25 = sum(abs(hs25).* ts25)/sum(abs(hs25));
rms_del_spr_25 = sqrt(sum(abs(hs25).*(ts25 - mean_del_25).^2)/sum(abs(hs25)));
figure(3)
plot(ts25(1:np25(1),1), 10*log10(abs(hs25(1:np25(1),1))))
hold on
p0_20 = 10*log10(PDP_0(:,10)/max(PDP_0(:, 10)));
pn5_n25 = 10*log10(PDP_n5(:,1)/max(PDP_n5(:,1)));
pp5_p25 = 10*log10(PDP_p5(:,11)/max(PDP_p5(:,11)));
plot(tau, p0_20)
plot(tau, pn5_n25)
plot(tau, pp5_p25)
hold off
legend('Simulated', '0 +20', '-5 -25', '+5 +25')
title('misalignment between 10-25 degrees')