% ************************************************************************
% This function reads the info about simulations from the file
% "par.inp" to get the simuls parameters

% Marion Thomas, last modified June 2017

% CALLS: 

%% ************************************************************************
function [fs,omega,beta,simuT] = getsimulparam(datadir)

fid = fopen([datadir 'Par.inp'], 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

%% get the values

%find value of fs
line = find(ismember(C{1}(:),'#---- Material parameters -------------- '));
str = char(C{1}(line+2));
ind1 = strfind(str,'fs=');
indtemp = strfind(str,',');
ind2 = indtemp(find(indtemp-ind1>0,1));
fschar=str(ind1+3:ind2-1);
indD=strfind(fschar,'d');
fschar(indD)='E';
fs = str2double(fschar); display(['fs = ', num2str(fs)]);

%find value of omega
str = char(C{1}(line+4));
ind1 = strfind(str,'omega=');
indtemp = strfind(str,',');
ind2 = indtemp(find(indtemp-ind1>0,1));
omegachar=str(ind1+6:ind2-1);
indD=strfind(omegachar,'d');
omegachar(indD)='E';
omega = str2double(omegachar); display(['omega = ', num2str(omega)]);

%find value of beta
str = char(C{1}(line+3));
ind1 = strfind(str,'beta=');
indtemp = strfind(str,'=');
ind2 = indtemp(find(indtemp-ind1>0,1));
betachar=str(ind2+1:end);
indD=strfind(betachar,'d');
betachar(indD)='E';
beta = str2double(betachar); display(['beta = ', num2str(beta)]);

%simulation duration
line = find(ismember(C{1}(:),'#---- Time scheme settings ---------------------- '));
str = char(C{1}(line+1));
ind1 = strfind(str,'TotalTime=');
ind2 = strfind(str,', courant =');
fschar=str(ind1+10:ind2-1);
simuT=str2double(fschar); 

