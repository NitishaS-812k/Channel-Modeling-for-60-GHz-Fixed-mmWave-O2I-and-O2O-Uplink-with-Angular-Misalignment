%%This script generates the pdp and TDL from Signal Analyzer data

%% On data : 190524-PHD_LAB-CESA-KONF1-CAL_SlotAnt.csv

%% Clearing all

%clear all;
close all;
%% System settings
% System constants
PlotType = 'IS';  % 'IS' ImageSc
CCFit = 'CIR';      % type of fitting

%% Numerical constants
AZnum = 13;      % [30,35,0,+-5,+-10,+-15,+-20,+-25]
ELnum = 3;        % [+-5,0]
Blim = [56,64];  % lower and upper cutoff frequency
Slices = (1:39);  % AMB - user defined CTFs
Elevation = 3;    %AMB
IntRat = 8;  %Interpolation ratio = 2^IntRat

%% Read file of measured data
[filename, path] = uigetfile('*.csv','Enter file name',' ');
[Freq, magSdBV, EL, AZ, AZxx, ELxx] = getSparCSV(strcat(path, filename)); 
ELstep = EL(2) - EL(1);
AZstep = AZ(2) - AZ(1);
dimdBV = size(magSdBV);

%% Frequency span reduction

% lower coordinate for freq. Blim
if (Blim(1)>Freq(1))
    NpL = (Blim(1) - Freq(1)) * 10 + 1;
else
    NpL = 1;
end

% upper coordinate for freq. Blim
if (Blim(2)<Freq(end))
    NpH = numel(Freq) - (Freq(end) - Blim(2))*10;
else
    NpH = numel(Freq);
end

% new frequency and data vectors
FScale = 10^9;
Freq = Freq(NpL:NpH,1)*FScale;
magSdBV = magSdBV(NpL:NpH,:);
freqL = Freq(1);
freqU = Freq(end);
Fs = (freqU-freqL)/(numel(Freq)-1);   %Frequency step in GHz

% FREQUENCY DOMAIN PROCESSING
%% convert measured data matrix to required shape
magSdB = zeros(81,3,13);
SizeMtrx = [13 3];  %[AZ EL]
for i = 1:numel(Freq)
    magSdB(i,:,:) = reshape (magSdBV(i,:), SizeMtrx)'; %magSdB is not used elsewhere
end


Np = NpH - NpL + 1;

%% Hilbert Transform
% phase calculation from antenna distance
H_Rx = 13.5;
H_Tx = 1.6;
D = 107;
%dAnt = sqrt(D^2 + (H_Rx-H_Tx)^2);
dAnt = 107.66;
TauCIR = dAnt/3e8;
magSV = 10.^(magSdBV./20);
SizemagSdBV = size(magSdBV);

% HT over all frequencies for particular AZ,EL
HphaseSRawV = zeros(SizemagSdBV(1),SizemagSdBV(2));
HphaseSV = zeros(SizemagSdBV(1),SizemagSdBV(2));
HSCmpxV = zeros(SizemagSdBV(1),SizemagSdBV(2));
for i = 1:SizemagSdBV(2)
    HphaseSRawV(:,i) = hilbtran(log(magSV(:,i))); % phase with HT only
    HphaseSV(:,i) = HphaseSRawV(:,i) + 2*pi*(Freq - Freq(1))*TauCIR; %phase with HT + Phase Comp.
    HSCmpxV(:,i) = (magSV(:,i).*exp(1j*HphaseSV(:,i))); %Complex CTF
end


%% Plots for phase comparison
% freq_points = 56:0.1:64;
% subplot(2,2,1);
% plot(freq_points,HphaseSRawV(:,19));
% xlabel('Frequency');
% ylabel('Phase');
% title('Phase Recovery using Hilbert Transform');
% 
%subplot(2,2,2);
%plot(freq_points,angle(HSCmpxV(:,19)));
%xlabel('Frequency');
%ylabel('Phase');
% title('Phase Recovery using HT + Phase Compensation');
% 
% subplot(2,2,3);
% plot(freq_points,HphaseSRawV(:,28));
% xlabel('Frequency');
% ylabel('Phase');
% title('Phase Recovery using Hilbert Transform');
% 
% subplot(2,2,4);
% plot(freq_points,HphaseSV(:,28));
% xlabel('Frequency');
% ylabel('Phase');
% title('Phase Recovery using HT + Phase Compensation');


%% IFFT

win(:,1) = ones(Np,1);
win(:,2) = hann (Np, 'periodic');
win(:,3) = hamming(Np, 'periodic');
win(:,4) = blackman(Np, 'periodic');
win(:,5) = flattopwin(Np,'periodic');

wintype = 3; %hamming windowing selected

HTCmpxV = zeros(SizemagSdBV(1),SizemagSdBV(2));
HmagTV = zeros(SizemagSdBV(1),SizemagSdBV(2));
for i=1:39
    HTCmpxV(:,i) = ifft (HSCmpxV(:,i).*win(:,wintype), Np, 1);
    HmagTV(:,i) = abs(HTCmpxV(:,i));
end

HmagdBT = 20*log10(HmagTV);    %magnitude log in time domain


%% RMS Delay Spread

SVar = magSV(:, Slices);    

SVarSize = size(SVar);

Tstop = 70;
BW = freqU - freqL;                 % bandwidth
Tr = abs(1./BW);
time = 0:Tr:Tr*(Np-1);
PDPsel = HmagTV(:, Slices).^2;       % Power Delay profile        
PDPselPks = zeros(SVarSize);      % Prepare array for MPC PDP
AvgTau = zeros(1, SVarSize(2));
RMSDelSpr = zeros(1, SVarSize(2));


for i = 1:SVarSize(2)
    % Find peaks of PDP and their locations
    [PksPDPsel, PksLocs] = findpeaks(PDPsel(:, i));
    
    % Find start time
    [MaxPk, PosPk] = max(PksPDPsel);    % Find max peak and its position
    Tstart = PksLocs(PosPk);
    TimeEval = Tstart:Tstop;            % Tstart:Tstop or 1:SVarSize(1);
    Tau = time(TimeEval);               
    
    % Set peaks in MPC PDP
    PDPselPks(PksLocs, i) = PksPDPsel;
    PDPselPks([1:Tstart - 1, Tstop + 1:SVarSize(1)], i) = 0;  % set zeros outside TimeEval
    
    % RMS Delay Spread
    PDP = PDPselPks;                    % PDPselPks or PDPsel
    AvgTau(i) = sum(PDP(TimeEval, i).* Tau')/sum(PDP(TimeEval, i));
    RMSDelSpr(i) = sqrt(sum(PDP(TimeEval, i).*(Tau - AvgTau(i))'.^2)/sum(PDP(TimeEval, i)));
end
PDPseldB = 10*log10(PDPsel);
PDPselPksdB = 10*log10(PDPselPks);

% RMS Delay Spread matrix
 indx=1;
 for i=1:3
     for j=1:13
         RMSmat(i,j) = RMSDelSpr(indx);
         indx = indx + 1; 
    end
 end
 
%% Modified PDP Plots

figure
plot((TimeEval-Tstart)*Tr*1e9,10*log10(PDPsel(TimeEval,19)));
hold on
plot((TimeEval-Tstart)*Tr*1e9,10*log10(PDPsel(TimeEval,1)));
plot((TimeEval-Tstart)*Tr*1e9,10*log10(PDPsel(TimeEval,10)));
plot((TimeEval-Tstart)*Tr*1e9,10*log10(PDPsel(TimeEval,39)));
xlabel('Delay w.r.t LOS (ns)');
ylabel('Power (in dB)');
title('Modified PDP for various spatial points');
legend('Point 19','Point 1','Point 10','Point 39')
hold off;


%% TDL Model

PDP_avg= mean(PDPsel,2);
 TDL_norm = zeros(1,10);
 TDL_norm_dB = zeros(1,10);
 seq = 10;
 t_samp=1; %Sampling rate(1 sample means 1*(1/BW) = 0.125 ns or 8 Gbits/s)
 N=10; %No of taps
 sample_pts = Tstart:t_samp:Tstart+(N-1)*t_samp;
 
 tseq=39;
 %for i=1:39  % Will be used to generate individual TDLs for various angles

 
 TDL_norm = PDP_avg(sample_pts)./PDP_avg(sample_pts(1)); %For Average TDL
 TDL_norm_dB = 10*log10(TDL_norm);

 %TDL_norm(i,:) = PDPsel(sample_pts,i)./PDPsel(sample_pts(1),i); % Will be used to generate individual TDLs for various angles
 %end
 

 %{
 figure;
 stem((sample_pts-Tstart)*Tr*1e9,TDL_norm);
 xlabel('Delay (ns)');
 ylabel('Normalised Tap value');
 title('Average TDL');
%}
 
 PDP_samp = TDL_norm;
 save PDP_samp.mat PDP_samp
 
%% For BER plots run the program BER_monteCarlo.m