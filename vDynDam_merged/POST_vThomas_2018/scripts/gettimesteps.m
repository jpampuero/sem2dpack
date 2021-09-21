% ************************************************************************
% This function reads the info about simulations from the file
% "sem2dpack.sta" to get the time stepping.

% Marion Thomas, last modified January 2017

% CALLS: 

%% ************************************************************************
function [time,outtime,steps] = gettimesteps(datadir)

fid = fopen([datadir 'sem2dpack.sta'], 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

line = find(ismember(C{1}(:),'Snapshot at timestep = 0'));
for i=0:2
str = char(C{1}(line-8+i));
ind = strfind(str,'=');
steps(i+1) = str2double(strtrim(str(ind+1:end)));
end

line = find(ismember(C{1}(:),'S n a p s h o t   O u t p u t s'));
str = char(C{1}(line+4));
ind = strfind(str,'=');
steps(i+2) = str2double(strtrim(str(ind+1:end)));


%time steps
dt = steps(1); 
%Number of time steps
steps(2)=steps(2)+1;
N = steps(2);
%Final time
totaltime = steps(3);
%outputs increment
outputsteps = steps(4);
%time
time = 0:dt:totaltime;%linspace(0,totaltime,N);
%time of outputs
outtime = (0:outputsteps:N)*dt;

end