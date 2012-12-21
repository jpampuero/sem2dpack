% SEM2D_READ_FAULT reads fault outputs from SEM2DPACK
%
% SYNTAX	data = sem2d_read_fault(name)
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
%		If output on each side of the fault (osides=T):
%  		data.d1t	displacement on side 1, fault parallel component
%  		data.d2t	displacement on side 2, fault parallel component
%  		data.d1n	displacement on side 1, fault normal component
%  		data.d2n	displacement on side 2, fault normal component
%  		data.v1t	velocity on side 1, fault parallel component
%  		data.v2t	velocity on side 2, fault parallel component
%  		data.v1n	velocity on side 1, fault normal component
%  		data.v2n	velocity on side 2, fault normal component
%
% NOTE		Fault normal components on each side of the fault are exported only in P-SV
%
function data = sem2d_read_fault(name)

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
%raw = reshape(raw(2:data.nx+1,:),[data.nx ndat data.nt]);
raw = reshape(raw(2:data.nx+1,:),data.nx,ndat,[]);

% Reformat each field [nx,nt]
data.d  = squeeze(raw(:,1,:)); 
data.v  = squeeze(raw(:,2,:)); 
data.st = squeeze(raw(:,3,:)); 
data.sn = squeeze(raw(:,4,:)); 
data.mu  = squeeze(raw(:,5,:)); 
if ndat==5+4
  data.d1t  = squeeze(raw(:,6,:)); 
  data.d2t  = squeeze(raw(:,7,:)); 
  data.v1t  = squeeze(raw(:,8,:)); 
  data.v2t  = squeeze(raw(:,9,:)); 
elseif ndat==5+4*2
  data.d1t  = squeeze(raw(:,6,:)); 
  data.d1n  = squeeze(raw(:,7,:)); 
  data.d2t  = squeeze(raw(:,8,:)); 
  data.d2n  = squeeze(raw(:,9,:)); 
  data.v1t  = squeeze(raw(:,10,:)); 
  data.v1n  = squeeze(raw(:,11,:)); 
  data.v2t  = squeeze(raw(:,12,:)); 
  data.v2n  = squeeze(raw(:,13,:)); 
end
