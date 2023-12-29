%This script generates the PDPs and the angular PDPs from O2O data

[filename3, path3] = uigetfile('*.csv','Enter file name',' ');
[Freq3, magSdBV3, EL3, AZ3, AZxx3, ELxx3] = getSparCSV(strcat(path3, filename3)); 
ELstep3 = EL3(2) - EL3(1);
AZstep3  = AZ3(2) - AZ3(1);
dimdBV3 = size(magSdBV3);
Freqrange3 = [56, 64];

% lower coordinate for freq. Blim
if (Freqrange3(1)>Freq3(1))
    NL3 = (Freqrange3(1) - Freq3(1)) * 10 + 1;
else
    NL3 = 1;
end

% upper coordinate for freq. Blim
if (Freqrange3(2)<Freq3(end))
    NH3 = numel(Freq3) - (Freq3(end) - Freqrange3(2))*10;
else
    NH3 = numel(Freq3);
end

Freq3 = Freq3(NL3:NH3,1)*(10^9);
magSdBV3 = magSdBV3(NL3:NH3,:);
freqL3 = Freq3(1);
freqU3 = Freq3(end);
Fs3 = (freqU3-freqL3)/(numel(Freq3)-1);   %Frequency step in GHz

NP3 = NH3-NL3+1;

%Hilbert transform phase calculation
dAnt3 = 98.1;
tauCIR3 = dAnt3/3e8;
magSV3 = 10.^(magSdBV3./20);
SizemagSdBV3 = size(magSdBV3);

%HT over all frequencies
HphaseSRawV3 = zeros(SizemagSdBV3);
HphaseSV3 = zeros(SizemagSdBV3);
HSCmpxV3 = zeros(SizemagSdBV3);
for i = 1:SizemagSdBV3(2)
    HphaseSRawV3(:,i) = hilbtran(log(magSV3(:,i))); % phase with HT only
    HphaseSV3(:,i) = HphaseSRawV3(:,i) + 2*pi*(Freq3 - Freq3(1))*tauCIR3; %phase with HT + Phase Comp.
    HSCmpxV3(:,i) = (magSV3(:,i).*exp(1j*HphaseSV3(:,i))); %Complex CTF
end
win3(:,1) = ones(NP3, 1);
win3(:,2) = hann(NP3, 'periodic');
win3(:,3) = hamming(NP3, 'periodic');
win3(:,4) = blackman(NP3, 'periodic');
win3(:,5) = flattopwin(NP3,'periodic');
wintype3 = 4;
HTCmpxV3 = zeros(SizemagSdBV3);
HmagTV3 = zeros(SizemagSdBV3);
for i = 1:SizemagSdBV3(2)
    HTCmpxV3(:,i) = ifft (HSCmpxV3(:,i).*win3(:,wintype3), NP3, 1);
    HmagTV3(:,i) = abs(HTCmpxV3(:,i));
end
HmagdBT3 = 20*log10(HmagTV3);    %magnitude log in time domain

%PDPs
PADP3 = HTCmpxV3.*conj(HTCmpxV3);
PADPdB3 = 10*log10(PADP3);
BW3 = freqU3 - freqL3;
Tr3 = abs(1./BW3);
time3 = 0:Tr3:Tr3*(NP3-1);
%PDPsel3 = HmagTV3(:, Slices).^2;       % Power Delay profile        
PDPselPks3 = zeros(SizemagSdBV3);      % Prepare array for MPC PDP
AvgTau3 = zeros(1, SizemagSdBV3(2));
RMSDelSpr3 = zeros(1, SizemagSdBV3(2));

for i = 1:SizemagSdBV3(2)
    %peaks and position of peaks at each position
    [peaks3, peaks_pos3] = findpeaks(PADP3(:,i));
    
    %LOS peak
    [Maxpeak3, pospeak3] = max(peaks3);
    tstart3 = peaks_pos3(pospeak3);
    timeeval3 = tstart3:81;
    tau3 = time3(timeeval3);
    % RMS Delay Spread
    PDP3 = PADP3;                    % PDPselPks or PDPsel
    AvgTau3(i) = sum(PDP3(timeeval3, i).* tau3')/sum(PDP3(timeeval3, i));
    RMSDelSpr3(i) = sqrt(sum(PDP3(timeeval3, i).*(tau3 - AvgTau3(i))'.^2)/sum(PDP3(timeeval3, i)));
end

PDPselPksdB3 = 10*log10(PDPselPks3);
timeee3 = (timeeval3-tstart3)*Tr3*1e9;

PAP3 = (1/NP3)*sum(PADP3,1);
PAPdB3 = 10*log10(PAP3);
el866 = find(ELxx3 == 8.6600);
el03 = find(ELxx3 == 0);
eln866 = find(ELxx3 == -8.6600);
eln13 = find(ELxx3 == -13.0000);
eln433 = find(ELxx3 == -4.3300);
el433 = find(ELxx3 == 4.3300);
A866 = AZxx3(el866);
A03 = AZxx3(el03);
An866 = AZxx3(eln866);
An13 = AZxx3(eln13);
An433 = AZxx3(eln433);
A433 = AZxx3(el433);

%PAP plot
figure(1)
plot(A866, PAPdB3(el866));
hold on
plot(A03, PAPdB3(el03));
plot(An866, PAPdB3(eln866));
plot(An13, PAPdB3(eln13));
plot(A433, PAPdB3(el433));
plot(An433, PAPdB3(eln433));
hold off
legend('Elevation +8.66', 'Elevation 0', 'Elevation -8.66', 'Elevation -13', 'Elevation +4.33', 'Elevation -4.33');
xlabel('azimuth angle');
ylabel('Power in dB');
title('Power Angle Profile')

%PDP plots
PDP_866 = PADP3(timeeval3, el866);     %Powers for elevation +8.66
PDP_866dB = PADPdB3(timeeval3, el866);
PDP_03 = PADP3(timeeval3, el03);   %Powers for elevation 0
PDP_03dB = PADPdB3(timeeval3, el03);
PDP_n866 = PADP3(timeeval3, eln866);   %Powers for elevation -8.66
PDP_n866dB = PADPdB3(timeeval3, eln866);
PDP_13 = PADP3(timeeval3, eln13);   %Powers for elevation -13
PDP_n13dB = PADPdB3(timeeval3, eln13);
PDP_433 = PADP3(timeeval3, el433);
PDP_433dB = PADPdB3(timeeval3, el433);
PDP_n433 = PADP3(timeeval3, eln433);
PDP_n433dB = PADPdB3(timeeval3, eln433);
%{
figure(2)
hold on
plot(timeee3, PDP_03dB(:,1)) %-25
plot(timeee3, PDP_03dB(:,11)) %+25
plot(timeee3, PDP_03dB(:,2)) %-20
plot(timeee3, PDP_03dB(:,10)) %+20
hold off
legend('-25', '+25', '-20', '+20')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation 0')

figure(3)
hold on
plot(timeee3, PDP_03dB(:,3)) %-15
plot(timeee3, PDP_03dB(:,9)) %+15
plot(timeee3, PDP_03dB(:,4)) %-10
plot(timeee3, PDP_03dB(:,8)) %+10
plot(timeee3, PDP_03dB(:,5)) %-5
plot(timeee3, PDP_03dB(:,7)) %+5
hold off
legend('-15', '+15', '-10', '+10', '-5', '+5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation 0')

figure(4)
hold on
plot(timeee3, PDP_n433dB(:,1)) %-25
plot(timeee3, PDP_n433dB(:,10)) %+25
plot(timeee3, PDP_n433dB(:,2)) %-20
plot(timeee3, PDP_n433dB(:,9)) %+20
hold off
legend('22.5', '-22.5', '17.5', '-17.5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation -4.33')

figure(5)
hold on
plot(timeee3, PDP_n433dB(:,3)) %-15
plot(timeee3, PDP_n433dB(:,8)) %+15
plot(timeee3, PDP_n433dB(:,4)) %-10
plot(timeee3, PDP_n433dB(:,7)) %+10
plot(timeee3, PDP_n433dB(:,5)) %-5
plot(timeee3, PDP_n433dB(:,6)) %+5
hold off
legend('12.5', '-12.5', '7.5', '-7.5', '2.5', '-2.5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation -4.33')

figure(6)
hold on
plot(timeee3, PDP_433dB(:,1)) %-25
plot(timeee3, PDP_433dB(:,10)) %+25
plot(timeee3, PDP_433dB(:,2)) %-20
plot(timeee3, PDP_433dB(:,9)) %+20
hold off
legend('22.5', '-22.5', '17.5', '-17.5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation +4.33')

figure(7)
hold on
plot(timeee3, PDP_433dB(:,3)) %-15
plot(timeee3, PDP_433dB(:,8)) %+15
plot(timeee3, PDP_433dB(:,4)) %-10
plot(timeee3, PDP_433dB(:,7)) %+10
plot(timeee3, PDP_433dB(:,5)) %-5
plot(timeee3, PDP_433dB(:,6)) %+5
hold off
legend('12.5', '-12.5', '7.5', '-7.5', '2.5', '-2.5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation +4.33')

figure(8)
hold on
plot(timeee3, PDP_n13dB(:,1)) %-25
plot(timeee3, PDP_n13dB(:,10)) %+25
plot(timeee3, PDP_n13dB(:,2)) %-20
plot(timeee3, PDP_n13dB(:,9)) %+20
hold off
legend('22.5', '-22.5', '17.5', '-17.5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation -13')

figure(9)
hold on
plot(timeee3, PDP_n13dB(:,3)) %-15
plot(timeee3, PDP_n13dB(:,8)) %+15
plot(timeee3, PDP_n13dB(:,4)) %-10
plot(timeee3, PDP_n13dB(:,7)) %+10
plot(timeee3, PDP_n13dB(:,5)) %-5
plot(timeee3, PDP_n13dB(:,6)) %+5
hold off
legend('12.5', '-12.5', '7.5', '-7.5', '2.5', '-2.5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation -13')

figure(10)
hold on
plot(timeee3, PDP_866dB(:,1)) %-25
plot(timeee3, PDP_866dB(:,11)) %+25
plot(timeee3, PDP_866dB(:,2)) %-20
plot(timeee3, PDP_866dB(:,10)) %+20
hold off
legend('-25', '+25', '-20', '+20')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation +8.66')

figure(11)
hold on
plot(timeee3, PDP_866dB(:,3)) %-15
plot(timeee3, PDP_866dB(:,9)) %+15
plot(timeee3, PDP_866dB(:,4)) %-10
plot(timeee3, PDP_866dB(:,8)) %+10
plot(timeee3, PDP_866dB(:,5)) %-5
plot(timeee3, PDP_866dB(:,7)) %+5
hold off
legend('-15', '+15', '-10', '+10', '-5', '+5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation +8.66')

figure(12)
hold on
plot(timeee3, PDP_n866dB(:,1)) %-25
plot(timeee3, PDP_n866dB(:,11)) %+25
plot(timeee3, PDP_n866dB(:,2)) %-20
plot(timeee3, PDP_n866dB(:,10)) %+20
hold off
legend('-25', '+25', '-20', '+20')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation -8.66')

figure(13)
hold on
plot(timeee3, PDP_n866dB(:,3)) %-15
plot(timeee3, PDP_n866dB(:,9)) %+15
plot(timeee3, PDP_n866dB(:,4)) %-10
plot(timeee3, PDP_n866dB(:,8)) %+10
plot(timeee3, PDP_n866dB(:,5)) %-5
plot(timeee3, PDP_n866dB(:,7)) %+5
hold off
legend('-15', '+15', '-10', '+10', '-5', '+5')
xlabel('delay in ns')
ylabel('power in dB')
title('PDPs for azimuths for elevation -8.66')
%}
%Angular spread values
TA03 = sum(A03.*PAP3(el03))/sum(PAP3(el03));
AS03 = sqrt(sum(((A03-TA03).^2).*(PAP3(el03)))/sum(PAP3(el03)));
TA866 = sum(A866.*PAP3(el866))/sum(PAP3(el866));
AS866 = sqrt(sum(((A866-TA866).^2).*(PAP3(el866)))/sum(PAP3(el866)));
TAn866 = sum(An866.*PAP3(eln866))/sum(PAP3(eln866));
ASn866 = sqrt(sum(((An866-TAn866).^2).*(PAP3(eln866)))/sum(PAP3(eln866)));
TAn13 = sum(An13.*PAP3(eln13))/sum(PAP3(eln13));
ASn13 = sqrt(sum(((An13-TAn13).^2).*(PAP3(eln13)))/sum(PAP3(eln13)));
TAn433 = sum(An433.*PAP3(eln433))/sum(PAP3(eln433));
ASn433 = sqrt(sum(((An433-TAn433).^2).*(PAP3(eln433)))/sum(PAP3(eln433)));
TA433 = sum(A433.*PAP3(el433))/sum(PAP3(el433));
AS433 = sqrt(sum(((A433-TA433).^2).*(PAP3(el433)))/sum(PAP3(el433)));
disp('Angular spread at elevation 0'); AS03
disp('Angular spread at elevation -8.66'); ASn866
disp('Angular spread at elevation +8.66'); AS866
disp('Angular spread at elevation -13'); ASn13
disp('Angular spread at elevation -4.33'); ASn433
disp('Angular spread at elevation +4.33'); AS433

figure(14)
scatter(A03, RMSDelSpr3(el03)*(1e9), '*')
hold on
scatter(A866, RMSDelSpr3(el866)*(1e9), '*')
scatter(An866, RMSDelSpr3(eln866)*(1e9), '*')
scatter(An13, RMSDelSpr3(eln13)*(1e9), '*')
scatter(An433, RMSDelSpr3(eln433)*(1e9),'*')
scatter(A433, RMSDelSpr3(el433)*(1e9), '*')
hold off
legend('Elevation 0', 'Elevation +8.64', 'Elevation -8.66', 'Elevation -13','Elevation -4.33', 'Elevation +4.33')
xlabel('Azimuth angle')
ylabel('RMS delay spread in ns')
title('RMS delay spread vs azimuth')

figure(15)
scatter(A03, AvgTau3(el03)*(1e9), '*')
hold on
scatter(A866, AvgTau3(el866)*(1e9), '*')
scatter(An866, AvgTau3(eln866)*(1e9), '*')
scatter(An13, AvgTau3(eln13)*(1e9), '*')
hold off
legend('Elevation 0', 'Elevation +8.64', 'Elevation -8.69', 'Elevation -17.36')
xlabel('Azimuth angle')
ylabel('RMS delay spread in ns')
title('Mean delay versus azimuth')

%{
figure(16)
[x11, y11] = idk(PDP_866dB, 5);
scatter(x11, y11)
title('-25(8.66)')
figure(17)
[x22, y22] = idk(PDP_866dB, 7);
scatter(x22, y22)
title('+25(8.66)')
figure(18)
[x33, y33] = idk(PDP_n866dB, 5);
scatter(x33, y33)
title('-25(-8.66)')
figure(19)
[x44, y44] = idk(PDP_n866dB, 7);
scatter(x44, y44)
title('25(-8.66)')
%}
figure(16)
%modified PAP
plot(A866, PAPdB3(el866));
hold on
plot(An866, PAPdB3(eln866));
plot(A433, PAPdB3(el433));
plot(An433, PAPdB3(eln433));
hold off
legend('Elevation +8.66', 'Elevation -8.66', 'Elevation +4.33', 'Elevation -4.33');
xlabel('azimuth angle');
ylabel('Power in dB');
title('Power Angle Profile')