% data = LoadFault(name)
%
% INPUT		name	[Flt01] prefix of header and data files (name_sem2d.*) 
%
% OUTPUT	data.nx		number of fault nodes
%		data.nt		number of time samples
%		data.dt		time step
%		data.x,data.z	coordinates of fault nodes
%		data.d		slip [nx,nt]
%		data.v		slip rate
%		data.st		shear stress
%		data.sn		normal stress
%		data.mu		friction coefficient
%
function data = LoadFault(name)

if ~exist('name','var'), name = 'Flt01'; end

% Read parameters from header file
hdr = strcat(name,'_sem2d.hdr');
if ~exist(hdr,'file'), data=[]; return, end
[data.nx,ndat,data.nt,data.dt] = textread(hdr,'%n%n%n%n',1,'headerlines',1);
[data.x,data.z] = textread(hdr,'%f%f','headerlines',4);

% Read fault data in a big matrix
dat  = strcat(name,'_sem2d.dat');
fid=fopen(dat); 
raw = fread(fid,[data.nx+2,inf],'single') ; 
fclose(fid);
raw = reshape(raw(2:data.nx+1,:),[data.nx ndat data.nt]);

% Reformat each field [nx,nt]
data.d  = squeeze(raw(:,1,:)); 
data.v  = squeeze(raw(:,2,:)); 
data.st = squeeze(raw(:,3,:)); 
data.sn = squeeze(raw(:,4,:)); 
data.mu  = squeeze(raw(:,5,:)); 
