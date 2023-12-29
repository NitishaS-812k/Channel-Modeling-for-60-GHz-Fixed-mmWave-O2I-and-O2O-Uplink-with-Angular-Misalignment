function [Freq, CIRData, EL, AZ, AZf, ELf] = getSparCSV(Filenum)
%% This function reformats the file 'Filenum' to return frequency vector in Freq, measured
% data in CIRData. All unique Elevation and Azimuthal angles are saved as
% EL, AZ respectively

%%
% p1 real part of s pram
% p2 Imaginary part of s pram


%% Initialize variables.
Delimiter = ';';
HeaderLines = 3;

%% Read data according to format string.

fileID = fopen(Filenum);
C = textscan(fileID,'%s');
fclose(fileID);

% generate vector of elevations
C1 = C{1};
ELstr = char(C1(2));
ELind = strfind(ELstr, ';');    % indicies of ';' position
ELstr(ELind) = ' ';             % ';' replacment by space
ELstr = ELstr(ELind(1)+1:end);  % first word removal 
ELf = sscanf(ELstr,'%f')';
EL = unique(sscanf(ELstr,'%f'), 'stable')';


% generate vector of azimuths
AZstr = char(C1(4));
AZind = strfind(AZstr, ';');    % indicies of ';' position
AZstr(AZind) = ' ';             % ';' replacment by space
AZstr = AZstr(AZind(1)+1:end);  % first word removal 
AZf = sscanf(AZstr,'%f')';
AZ = unique(sscanf(AZstr,'%f'), 'stable')';

dataArray = dlmread(Filenum, Delimiter, HeaderLines, 0);


%% Close the text file.

%% Allocate imported array to column variable names
Freq = dataArray(:, 1);
CIRData = dataArray(:, 2:end);


%% Clear temporary variables
clearvars Delimiter HeaderLines dataArray;