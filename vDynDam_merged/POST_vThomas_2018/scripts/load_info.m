% Script to download information from the output folder

% Marion Thomas last modified July 2018

%CALLS: 
clc
%==========================================================================
%% 
disp(' '); 
disp('========================================');
disp('Bulk Properties')
disp('========================================');

%Elastic Bulk properties
disp(' '); disp('%%%%Elastic model');
[infoE,NmatE] = getsimulinfoE(datadirE);
disp(' '); disp('%%%%Flat Elastic model');
[infoF,NmatF] = getsimulinfoE(datadirF);


%Damage Bulk properties
disp(' '); disp('%%%%Damage model');
[fric,omega,beta,simuD] = getsimulparam(datadirD);
[simuE] = getsimulparamE(datadirE);
[simuF] = getsimulparamE(datadirF);
[infoD,MuS,MuD,ndof,NmatD] = getsimulinfo(datadirD);

%Time
[tallD,outtimeD,stepsD] = gettimesteps(datadirD);
dtD=stepsD(1);
TD=outtimeD;
[tallE,outtimeE,stepsE] = gettimesteps(datadirE);
dtE=stepsE(1);
TE=outtimeE;
[tallF,outtimeF,stepsF] = gettimesteps(datadirF);
dtF=stepsF(1);
TF=outtimeF;

%Snapshot
ooD=000:001:floor(stepsD(2)/stepsD(4));
ooE=000:001:floor(stepsE(2)/stepsE(4));
ooF=000:001:floor(stepsF(2)/stepsF(4));

%Info for plots
if NmatD==1; infoP=[infoD(1:4),infoD(1:4)]; else infoP=[infoD(1,1:4),infoD(end,1:4)];end
dt = (tplot-tstart)/ndt; %time steps for the slip and slip velocity plots



