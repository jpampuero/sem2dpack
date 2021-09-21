% ************************************************************************
% This function reads the info about simulations from the file
% "par.inp" to get the simuls parameters

% Marion Thomas, last modified June 2017

% CALLS: 

%% ************************************************************************
function [simuT] = getsimulparamE(datadir)

fid = fopen([datadir 'Par.inp'], 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

%% get the values

%simulation duration
line = find(ismember(C{1}(:),'#---- Time scheme settings ---------------------- '));
str = char(C{1}(line+1));
ind1 = strfind(str,'TotalTime=');
ind2 = strfind(str,', courant =');
fschar=str(ind1+10:ind2-1);
simuT=str2double(fschar); 

