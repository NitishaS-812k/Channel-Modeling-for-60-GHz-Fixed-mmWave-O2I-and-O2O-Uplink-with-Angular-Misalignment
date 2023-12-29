%Generates the PADP for O2I
PADP = HTCmpxV.*conj(HTCmpxV);
k = size(PADP);
nfreq = k(1);    %number of frequency points
npos = k(2);     %number of positions
PAP = (1/nfreq)*sum(PADP, 1);  %power angle profile
PAPdB = 10*log10(PAP);
el0 = find(ELxx == 0);   %indices with elevation = 0
eln5 = find(ELxx == -5); %indices with elevation = -5
elp5 = find(ELxx == 5);  %indices with elevation = +5
A0 = AZxx(el0);
An5 = AZxx(eln5);
Ap5 = AZxx(elp5);
%power angle profile plot for different elevations
figure(1)
plot(A0, PAPdB(el0));
hold on
plot(An5, PAPdB(eln5));
plot(Ap5, PAPdB(elp5));
hold off
legend('Elevation 0', 'Elevation -5', 'Elevation +5');
xlabel('azimuth angle');
ylabel('Power in dB');
title('Power Angle Profile')
Maxpeaks = zeros(1, npos);
peaks_all = zeros(npos, 61);
for i = 1:npos
    %peaks and position of peaks at each position
    [peaks, peaks_pos] = findpeaks(PADP(:,i));
    peaks_all(i, peaks_pos) = 10*log10(peaks);
    %LOS peak
    [Maxpeak, pospeak] = max(peaks);
    Maxpeaks(i) = Maxpeak;
    tstart = peaks_pos(pospeak);
    timeeval = tstart:70;
end
%m = max(PADP(timeeval, 19));
%PADPnorm = PADP./m;
PADPdB = 10*log10(PADP);
tau = (timeeval-tstart)*Tr*1e9;
PDP_0 = PADP(timeeval, el0);     %Powers for elevation 0
PDP_0dB = PADPdB(timeeval, el0);
PDP_n5 = PADP(timeeval, eln5);   %Powers for elevation -5
PDP_n5dB = PADPdB(timeeval, elp5);
PDP_p5 = PADP(timeeval, elp5);   %Powers for elevation +5
PDP_p5dB = PADPdB(timeeval, elp5);
%{
figure(2)
hold on
plot(tau, PDP_0dB(:, 1)) %-25
plot(tau, PDP_0dB(:, 11))
plot(tau, PDP_0dB(:, 2)) %-20
plot(tau, PDP_0dB(:, 10))
plot(tau, PDP_0dB(:, 3)) %-15
plot(tau, PDP_0dB(:, 9))
hold off
legend('-25','25','-20','20','-15','15')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation 0')

figure(3)
hold on
plot(tau, PDP_0dB(:, 4)) %-10
plot(tau, PDP_0dB(:, 8))
plot(tau, PDP_0dB(:, 5)) %-5
plot(tau, PDP_0dB(:, 7))
hold off
legend('-10', '+10', '-5', '+5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation 0')

figure(4)
hold on
plot(tau, PDP_n5dB(:, 1)) %-25
plot(tau, PDP_n5dB(:, 11))
plot(tau, PDP_n5dB(:, 2)) %-20
plot(tau, PDP_n5dB(:, 10))
plot(tau, PDP_n5dB(:, 3)) %-15
plot(tau, PDP_n5dB(:, 9))
hold off
legend('-25','25','-20','20','-15','15')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation -5')

figure(5)
hold on
plot(tau, PDP_n5dB(:, 4)) %-10
plot(tau, PDP_n5dB(:, 8))
plot(tau, PDP_n5dB(:, 5)) %-5
plot(tau, PDP_n5dB(:, 7))
plot(tau, PDP_n5dB(:, 6))
hold off
legend('-10', '+10', '-5', '+5', '0')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation -5')

figure(6)
hold on
plot(tau, PDP_p5dB(:, 1)) %-25
plot(tau, PDP_p5dB(:, 11))
plot(tau, PDP_p5dB(:, 2)) %-20
plot(tau, PDP_p5dB(:, 10))
plot(tau, PDP_p5dB(:, 3)) %-15
plot(tau, PDP_p5dB(:, 9))
hold off
legend('-25','25','-20','20','-15','15')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation +5')

figure(7)
hold on
plot(tau, PDP_n5dB(:, 4)) %-10
plot(tau, PDP_n5dB(:, 8))
plot(tau, PDP_n5dB(:, 5)) %-5
plot(tau, PDP_n5dB(:, 7))
plot(tau, PDP_n5dB(:, 6))
hold off
legend('-10', '+10', '-5', '+5', '0')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for positive azimuths for elevation +5')

az30 = find(AZxx == 30);
PDP_az30 = PADP(timeeval, az30);
az25 = find(AZxx == 25);
PDP_az25 = PADP(timeeval, az25);
az35 = find(AZxx == 35);
PDP_az35 = PADP(timeeval, az35);
az20 = find(AZxx == 20);
PDP_az20 = PADP(timeeval, az20);
az15 = find(AZxx == 15);
PDP_az15 = PADP(timeeval, az15);
az10 = find(AZxx == 10);
PDP_az10 = PADP(timeeval, az10);
az5 = find(AZxx == 5);
PDP_az5 = PADP(timeeval, az5);
az0 = find(AZxx == 0);
PDP_az0 = PADP(timeeval, az0);
azn5 = find(AZxx == -5);
PDP_azn5 = PADP(timeeval, azn5);
azn10 = find(AZxx == -10);
PDP_azn10 = PADP(timeeval, azn10);
azn15 = find(AZxx == -15);
PDP_azn15 = PADP(timeeval, azn15);
azn20 = find(AZxx == -20);
PDP_azn20 = PADP(timeeval, azn20);
azn25 = find(AZxx == -25);
PDP_azn25 = PADP(timeeval, azn25);

figure(8)
hold on
for ii = 1:length(az30)
    plot(tau, 10*log10(PDP_az30(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(az30(1))), num2str(ELxx(az30(2))), num2str(ELxx(az30(3))))
title('PDPs for azimuth 30')

figure(9)
hold on
for ii = 1:length(az20)
    plot(tau, 10*log10(PDP_az20(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(az20(1))), num2str(ELxx(az20(2))), num2str(ELxx(az20(3))))
title('PDPs for azimuth 20')

figure(10)
hold on
for ii = 1:length(az25)
    plot(tau, 10*log10(PDP_az25(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(az25(1))), num2str(ELxx(az25(2))), num2str(ELxx(az25(3))))
title('PDPs for azimuth 25')

figure(11)
hold on
for ii = 1:length(az35)
    plot(tau, 10*log10(PDP_az35(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(az35(1))), num2str(ELxx(az35(2))), num2str(ELxx(az35(3))))
title('PDPs for azimuth 35')

figure(12)
hold on
for ii = 1:length(az15)
    plot(tau, 10*log10(PDP_az15(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(az15(1))), num2str(ELxx(az15(2))), num2str(ELxx(az15(3))))
title('PDPs for azimuth 15')

figure(13)
hold on
for ii = 1:length(az10)
    plot(tau, 10*log10(PDP_az10(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(az10(1))), num2str(ELxx(az10(2))), num2str(ELxx(az10(3))))
title('PDPs for azimuth 10')

figure(14)
hold on
for ii = 1:length(az5)
    plot(tau, 10*log10(PDP_az5(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(az5(1))), num2str(ELxx(az5(2))), num2str(ELxx(az5(3))))
title('PDPs for azimuth 5')

figure(15)
hold on
for ii = 1:length(az0)
    plot(tau, 10*log10(PDP_az0(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(az0(1))), num2str(ELxx(az0(2))), num2str(ELxx(az0(3))))
title('PDPs for azimuth 0')

figure(16)
hold on
for ii = 1:length(azn5)
    plot(tau, 10*log10(PDP_azn5(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(azn5(1))), num2str(ELxx(azn5(2))), num2str(ELxx(azn5(3))))
title('PDPs for azimuth -5')

figure(17)
hold on
for ii = 1:length(azn10)
    plot(tau, 10*log10(PDP_azn10(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(azn10(1))), num2str(ELxx(azn10(2))), num2str(ELxx(azn10(3))))
title('PDPs for azimuth -10')

figure(18)
hold on
for ii = 1:length(azn15)
    plot(tau, 10*log10(PDP_azn15(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(azn15(1))), num2str(ELxx(azn15(2))), num2str(ELxx(azn15(3))))
title('PDPs for azimuth -15')

figure(19)
hold on
for ii = 1:length(az10)
    plot(tau, 10*log10(PDP_azn20(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(azn20(1))), num2str(ELxx(azn20(2))), num2str(ELxx(azn20(3))))
title('PDPs for azimuth -20')

figure(20)
hold on
for ii = 1:length(azn25)
    plot(tau, 10*log10(PDP_azn25(:,ii)))
end
hold off
xlabel('Delay in ns')
ylabel('Power in dB')
legend(num2str(ELxx(azn25(1))), num2str(ELxx(azn25(2))), num2str(ELxx(azn25(3))))
title('PDPs for azimuth -25')

PAP_az30 = PAPdB(az30);
PAP_az0 = PAPdB(az0);
PAP_az35 = PAPdB(az35);
PAP_az25 = PAPdB(az25);
PAP_az20 = PAPdB(az20);
PAP_az15 = PAPdB(az15);
PAP_az10 = PAPdB(az10);
PAP_az5 = PAPdB(az5);
PAP_azn10 = PAPdB(azn10);
PAP_azn15 = PAPdB(azn15);
PAP_azn20 = PAPdB(azn20);
PAP_azn25 = PAPdB(azn25);
PAP_azn5 = PAPdB(azn5);
figure(21)
hold on
plot(ELxx(az35), PAP_az35)
plot(ELxx(az30),PAP_az30)
plot(ELxx(az25), PAP_az25)
plot(ELxx(az20), PAP_az20)
plot(ELxx(az15), PAP_az15)
plot(ELxx(az10), PAP_az10)
plot(ELxx(az5), PAP_az5)
plot(ELxx(az0), PAP_az0)
plot(ELxx(azn5), PAP_azn5)
plot(ELxx(azn10), PAP_azn10)
plot(ELxx(azn15), PAP_azn15)
plot(ELxx(azn20), PAP_azn20)
plot(ELxx(azn25), PAP_azn25)
hold off
legend('35','30','25','20','15','10','5','0','-5','-10','-15','-20','-25')
xlabel('Elevation angle')
ylabel('Power in dB')
title('Elevation power profile')

%{
figure(22)
[x1, y1] = idk(PDP_0dB, 5); %-25
scatter(x1, y1)
title('-5_0')
figure(23)
[x2, y2] = idk(PDP_0dB, 7); %+25
scatter(x2, y2)
title('+5_0')
figure(24)
[x3, y3] = idk(PDP_n5dB, 5); %-10
scatter(x3, y3)
title('-5_n5')
figure(25)
[x4, y4] = idk(PDP_n5dB, 7); %+10
scatter(x4, y4)
title('+5_n5')
figure(26)
[x5,y5] = idk(PDP_p5dB, 5); %-20
scatter(x5, y5)
title('-5_p5')
figure(27)
[x6, y6] = idk(PDP_p5dB, 7); %+20
scatter(x6, y6)
title('+5_p5')
%}
%}
figure(22)
%modified PAP plots
hold on
plot(An5, PAPdB(eln5));
plot(Ap5, PAPdB(elp5));
hold off
legend('Elevation -5', 'Elevation +5');
xlabel('azimuth angle');
ylabel('Power in dB');
title('Power Angle Profile')