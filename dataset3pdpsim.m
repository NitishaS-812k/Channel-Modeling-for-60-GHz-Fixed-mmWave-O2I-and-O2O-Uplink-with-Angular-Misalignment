%SV Parameters gamma(g), Gamma(G), lambda(l) and Lambda(L) sorted by
%elevation(0, 8.66, -8.66, 4.33, -4.33, -13)
g1_0 = cl13(el03);     g2_0 = cl23(el03);   l1_0 = l13(el03);        l2_0 = l23(el03);
g1_866 = cl13(el866);  g2_866 = cl23(el866); l1_866 = l13(el866);    l2_866 = l23(el866);
g1_n866 = cl13(eln866); g2_n866 = cl23(eln866); l1_n866 = l13(eln866);l2_n866 = l23(eln866);
g1_n13 = cl13(eln13);g2_n13 = cl23(eln13); l1_n13 = l13(eln13);l2_n13= l23(eln13);
g1_433 = cl13(el433); g2_433 = cl23(el433); l1_433 = l13(el433); l2_433 = l23(el433);
g1_n433 = cl13(eln433); g2_n433 = cl23(eln433); l1_n433 = l13(eln433); l2_n433 = l23(eln433);
G1_0 = GAMMA3(el03);     L1_0 = LAMBDA3(el03);  g3_0 = cl33(el03); l3_0 = l33(el03);
G1_866 = GAMMA3(el866);  L1_866 = LAMBDA3(el866); g3_866 = cl33(el866); l3_866 = l33(el866);
G1_n866 = GAMMA3(eln866);  L1_n866 = LAMBDA3(eln866); g3_n866 = cl33(eln866); l3_n866 = l33(eln866);
G1_n13 = GAMMA3(eln13);    L1_n13 = LAMBDA3(eln13); g3_n13 = cl33(eln13); l3_n13 = l33(eln13);
G1_433 = GAMMA3(el433); L1_433 = LAMBDA3(el433); g3_433 = cl33(el433); l3_433 = l33(el433);
G1_n433 = GAMMA3(eln433);L1_n433 = LAMBDA3(eln433); g3_n433 = cl33(eln433); l3_n433 = l33(eln433);

%misalignment of 10
mis10_0 = find(A03 >= -10 & A03 <= 10);
los = mis10_0(3);
mis10_0(3) = [];
mis10_866 = find(A866 >= -10 & A866 <= 10); 
mis10_n866 = find(An866 >= -10 & An866 <= 10); 
mis10_n13 = find(An13 >= -10 & An13 <= 10);
mis10_433 = find(A433 >= -10 & A433 <= 10);
mis10_n433 = find(An433 >= -10 & An433 <= 10);

g1_m10 = 10^(mean([g1_0(mis10_0)', g1_866(mis10_866)', g1_n866(mis10_n866)', g1_n13(mis10_n13)', g1_433(mis10_433)', g1_n433(mis10_n433)'])/10);
%g2_m10 = 0.6;
%g3_m10 = 0.5;
 g2_m10 = 10^(mean([g2_0(mis10_0)', g2_866(mis10_866)', g2_n866(mis10_n866)', g2_n13(mis10_n13)', g2_433(mis10_433)', g2_n433(mis10_n433)'])/10);
 g3_m10 = 10^(mean([g3_0(mis10_0)', g3_866(mis10_866)', g3_n866(mis10_n866)', g3_n13(mis10_n13)', g3_433(mis10_433)', g3_n433(mis10_n433)'])/10);
l1_m10 = mean([l1_0(mis10_0)', l1_866(mis10_866)', l1_n866(mis10_n866)', l1_n13(mis10_n13)', l1_433(mis10_433)', l1_n433(mis10_n433)']);
l2_m10 = mean([l2_0(mis10_0)', l2_866(mis10_866)', l2_n866(mis10_n866)', l2_n13(mis10_n13)', l2_433(mis10_433)', l2_n433(mis10_n433)']);
l3_m10 = mean([l3_0(mis10_0)', l3_866(mis10_866)', l3_n866(mis10_n866)', l3_n13(mis10_n13)', l3_433(mis10_433)', l3_n433(mis10_n433)']);
 G1_m10 = 10^(mean([G1_0(mis10_0)', G1_866(mis10_866)', G1_n866(mis10_n866)', G1_n13(mis10_n13)', G1_433(mis10_433)', G1_n433(mis10_n433)'])/10);
L1_m10 = mean([L1_0(mis10_0)', L1_866(mis10_866)', L1_n866(mis10_n866)', L1_n13(mis10_n13)', L1_433(mis10_433)', L1_n433(mis10_n433)']);
%G1_m10 = 0.9;
%misalignment of 10 - 25 
mis25_0 = [find(A03 > 10 & A03 <= 25), find(A03 < -10 & A03 >= -25)];
mis25_866 = [find(A866 > 10 & A866 <= 25), find(A866 < -10 & A866 >= -25)];
mis25_n866 = [find(An866 > 10 & An866 <= 25), find(An866 < -10 & An866 >= -25)];
mis25_n13 = [find(An13 > 10 & An13 <= 25), find(An13 < -10 & An13 >= -25)];
mis25_433 = [find(A433 > 10 & A433 <= 25), find(A433 < -10 & A433 >= -25)];
mis25_n433 = [find(An433 > 10 & An433 <= 25), find(A433 < -10 & An433 >= -25)];

g1_m25 = 10^(mean([g1_0(mis25_0)', g1_866(mis25_866)', g1_n866(mis25_n866)', g1_n13(mis25_n13)', g1_433(mis25_433)', g1_n433(mis25_n433)'])/10);
%g2_m25 = 0.7;
%g3_m25 = 0.6;
g2_m25 = 10^(mean([g2_0(mis25_0)', g2_866(mis25_866)', g2_n866(mis25_n866)', g2_n13(mis25_n13)', g2_433(mis25_433)', g2_n433(mis25_n433)'])/10);
g3_m25 = 10^(mean([g3_0(mis25_0)', g3_866(mis25_866)', g3_n866(mis25_n866)', g3_n13(mis25_n13)', g3_433(mis25_433)', g3_n433(mis25_n433)'])/10);
l1_m25 = mean([l1_0(mis25_0)', l1_866(mis25_866)', l1_n866(mis25_n866)', l1_n13(mis25_n13)', l1_433(mis25_433)', l1_n433(mis25_n433)']);
l2_m25 = mean([l2_0(mis25_0)', l2_866(mis25_866)', l2_n866(mis25_n866)', l2_n13(mis25_n13)', l2_433(mis25_433)', l2_n433(mis25_n433)']);
l3_m25 = mean([l3_0(mis25_0)', l3_866(mis25_866)', l3_n866(mis25_n866)', l3_n13(mis25_n13)', l3_433(mis25_433)', l3_n433(mis25_n433)']);
G1_m25 = 10^(mean([G1_0(mis25_0)', G1_866(mis25_866)', G1_n866(mis25_n866)', G1_n13(mis25_n13)', G1_433(mis25_433)', G1_n433(mis25_n433)'])/10);
L1_m25 = mean([L1_0(mis25_0)', L1_866(mis25_866)', L1_n866(mis25_n866)', L1_n13(mis25_n13)', L1_433(mis25_433)', L1_n433(mis25_n433)']);
%G1_m25 = 0.99;
%LOS
los01 = 10^(g1_0(los)/10); %gamma1 for los
los02 = 10^(g2_0(los)/10); %gamma2 for los
los03 = 10^(g3_0(los)/10); %gamma3 for los
l_los1 =  l1_0(los);%lambda 1 for los
l_los2 = l2_0(los);%lambda 2 for los
l_los3 = l3_0(los); %lambda3 for los
G_los = 10^(G1_0(los)/10); 
L_los = L1_0(los);
[p00,t00,~,np00] = pdp_svsim_og(L_los,[l_los1, l_los2, l_los3],G_los,[los01,los02, los03],3,1,1,3,0,1);
figure(1)
plot(t00(1:np00(1),1), 10*log10(abs(p00(1:np00(1),1))))
hold on
p0 = 10*log10(PDP_03(:,6)/max(PDP_03(:, 6)));
plot(timeee3, p0)
hold off
legend('Simulated', '0_0')
title('LOS')

%test PDPs
%for PDPs in 0-10 misalignment range
[hs10,ts10,t0s10,np10] = pdp_svsim_og(L1_m10,[l1_m10, l2_m10, l3_m10],G1_m10,[g1_m10,g2_m10,g3_m10],3,1,1,3,0,1);
mean_del_10 = sum(abs(hs10).* ts10)/sum(abs(hs10));
rms_del_spr_10 = sqrt(sum(abs(hs10).*(ts10 - mean_del_10).^2)/sum(abs(hs10)));
figure(2)
plot(ts10(1:np10(1),1), 10*log10(abs(hs10(1:np10(1),1))))
hold on
p0_n5 = 10*log10(PDP_03(:,5)/max(PDP_03(:,5)));
p10_866 = 10*log10(PDP_866(:,8)/max(PDP_866(:,8)));
pn10_n866 = 10*log10(PDP_n866(:,4)/max(PDP_n866(:,4)));
p7_n433= 10*log10(PDP_n433(:,4)/max(PDP_n433(:,4)));
pn7_433 = 10*log10(PDP_433(:,7)/max(PDP_433(:,7)));
p2_n13 = 10*log10(PDP_13(:,5)/max(PDP_13(:,5)));
plot(timeee3, p0_n5)
plot(timeee3, p10_866)
plot(timeee3, pn10_n866)
plot(timeee3, p7_n433)
plot(timeee3, pn7_433)
plot(timeee3, p2_n13)
hold off
legend('Simulated','0 -5','+8.66 +10','-8.66 -10','-4.33 +7.5', '+4.33 -7.5', '-13 +2.5')
title('misalignment between 0-10 degrees')

%for PDPs in 10-25 misalignment range
[hs25,ts25,t0s25,np25] = pdp_svsim_og(L1_m25,[l1_m25,l2_m25,l3_m25],G1_m25,[g1_m25,g2_m25,g3_m25],3,1,1,3,0,1);
mean_del_25 = sum(abs(hs25).* ts25)/sum(abs(hs25));
rms_del_spr_25 = sqrt(sum(abs(hs25).*(ts25 - mean_del_25).^2)/sum(abs(hs25)));
figure(3)
plot(ts25(1:np25(1),1), 10*log10(abs(hs25(1:np25(1),1))))
hold on
p0_20 = 10*log10(PDP_03(:,10)/max(PDP_03(:, 10)));
pn22_n13 = 10*log10(PDP_13(:,10)/max(PDP_13(:,10)));
p25_866 = 10*log10(PDP_866(:,11)/max(PDP_866(:,11)));
pn15_n866 = 10*log10(PDP_n866(:,3)/max(PDP_n866(:,3)));
p7_433 = 10*log10(PDP_433(:,4)/max(PDP_433(:,4)));
pn7_n433 = 10*log10(PDP_n433(:, 7)/max(PDP_n433(:,7)));
plot(timeee3, p0_20)
plot(timeee3, pn22_n13)
plot(timeee3, p25_866)
plot(timeee3, pn15_n866)
plot(timeee3, p7_433)
plot(timeee3, pn7_n433)
hold off
legend('Simulated', '0 +20', '-13 -22.5', '+8.66 +25', '-8.66 -15', '+4.33 +7.5', '-4.33 -7.5')
title('misalignment between 10-25 degrees')